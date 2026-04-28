//! .vpf (VariantPaths File) reader.  See build_vpf.py for the writer.

use anyhow::{anyhow, bail, Result};
use byteorder::{LittleEndian, ReadBytesExt};
use std::fs::File;
use std::io::{Cursor, Read};
use std::path::Path;

pub const MAGIC: &[u8; 4] = b"VPF1";
pub const HEADER_SZ: usize = 64;
pub const RECORD_SZ: usize = 32;
pub const FLAG_ZSTD: u16 = 1 << 0;
pub const NIL_IDX: u32 = 0xFFFF_FFFF;

#[derive(Debug, Clone)]
pub struct Header {
    pub version: u16,
    pub flags: u16,
    pub n_chroms: u32,
    pub n_samples: u32,
    pub n_classes: u32,
    pub n_bubbles: u64,
    pub built_unix: u64,
    pub reference_id: String,
}

#[derive(Debug, Clone)]
pub struct Pool {
    pub strings: Vec<String>,
}

impl Pool {
    fn read(rdr: &mut Cursor<&[u8]>) -> Result<Self> {
        let size = rdr.read_u32::<LittleEndian>()? as usize;
        let n = rdr.read_u32::<LittleEndian>()? as usize;
        let mut offsets = Vec::with_capacity(n);
        for _ in 0..n {
            offsets.push(rdr.read_u32::<LittleEndian>()? as usize);
        }
        let pos = rdr.position() as usize;
        let inner = rdr.get_ref();
        if pos + size > inner.len() {
            bail!("string pool size {} exceeds buffer", size);
        }
        let block = &inner[pos..pos + size];
        rdr.set_position((pos + size) as u64);

        let mut strings = Vec::with_capacity(n);
        for &off in &offsets {
            if off >= block.len() { bail!("pool offset {} out of range", off); }
            // NUL-terminated
            let end = block[off..].iter().position(|&b| b == 0).unwrap_or(block.len() - off);
            let s = std::str::from_utf8(&block[off..off + end])
                .map_err(|e| anyhow!("invalid utf8 in pool: {}", e))?;
            strings.push(s.to_string());
        }
        Ok(Pool { strings })
    }

    pub fn get(&self, idx: u32) -> Option<&str> {
        if idx == NIL_IDX { return None; }
        self.strings.get(idx as usize).map(|s| s.as_str())
    }
}

#[derive(Debug, Clone)]
pub struct ChromIdx {
    pub chrom_name_idx: u32,
    pub length_bp: u32,
    pub bubble_offset: u64,
    pub bubble_count: u32,
}

/// 32-byte fixed record.
#[derive(Debug, Clone, Copy)]
pub struct BubbleRec {
    pub start: i32,
    pub end: i32,
    pub vaf_q16: u16,
    pub chrom_idx: u16,
    pub length_log10_q: u8,
    pub flags: u8,
    pub sample_idx: u8,
    pub class_idx: u8,
    pub n_alts: u8,
    pub total_reads_log: u8,
    pub dbvar_recip_q8: u8,
    pub bubble_name_idx: u32,
    pub dbvar_id_idx: u32,
}

impl BubbleRec {
    pub fn vaf(&self) -> f32 { (self.vaf_q16 as f32) / 65535.0 }
    pub fn dbvar_recip(&self) -> f32 { (self.dbvar_recip_q8 as f32) / 255.0 }
    pub fn length(&self) -> i32 { (self.end - self.start).max(0) }
    /// approximate (since we quantized total_reads_log = log2(R)*8)
    pub fn total_reads_approx(&self) -> u32 {
        if self.total_reads_log == 0 { 0 }
        else { (2.0_f32).powf(self.total_reads_log as f32 / 8.0) as u32 }
    }
    /// `is_shared` = bubble belongs to a locus that recurs in ≥2 samples
    /// (orthogonal to the primary class).  Encoded in `flags` bit 0.
    pub fn is_shared(&self) -> bool { self.flags & 0x01 != 0 }
}

#[derive(Debug, Clone)]
pub struct Bbf {
    pub header: Header,
    pub samples: Pool,
    pub classes: Pool,
    pub chroms: Pool,
    pub bubble_names: Pool,
    pub dbvar_ids: Pool,
    pub chrom_index: Vec<ChromIdx>,
    pub bubbles: Vec<BubbleRec>,
}

impl Bbf {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut buf = Vec::new();
        File::open(path.as_ref())?.read_to_end(&mut buf)?;
        Self::parse(&buf)
    }

    pub fn parse(buf: &[u8]) -> Result<Self> {
        if buf.len() < HEADER_SZ { bail!("file too small for header"); }
        let mut rdr = Cursor::new(&buf[..HEADER_SZ]);
        let mut magic = [0u8; 4];
        rdr.read_exact(&mut magic)?;
        if &magic != MAGIC { bail!("bad magic: {:?}", magic); }
        let version = rdr.read_u16::<LittleEndian>()?;
        let flags = rdr.read_u16::<LittleEndian>()?;
        let n_chroms = rdr.read_u32::<LittleEndian>()?;
        let n_samples = rdr.read_u32::<LittleEndian>()?;
        let n_classes = rdr.read_u32::<LittleEndian>()?;
        let n_bubbles = rdr.read_u64::<LittleEndian>()?;
        let built_unix = rdr.read_u64::<LittleEndian>()?;
        let mut ref_bytes = [0u8; 16];
        rdr.read_exact(&mut ref_bytes)?;
        let mut _reserved = [0u8; 12];
        rdr.read_exact(&mut _reserved)?;
        let reference_id = std::str::from_utf8(&ref_bytes)
            .unwrap_or("")
            .trim_end_matches('\0')
            .to_string();

        let header = Header {
            version, flags,
            n_chroms, n_samples, n_classes, n_bubbles,
            built_unix, reference_id,
        };

        // Body
        let body = if flags & FLAG_ZSTD != 0 {
            zstd::decode_all(&buf[HEADER_SZ..])
                .map_err(|e| anyhow!("zstd decode failed: {}", e))?
        } else {
            buf[HEADER_SZ..].to_vec()
        };

        let mut br = Cursor::new(body.as_slice());

        let samples = Pool::read(&mut br)?;
        let classes = Pool::read(&mut br)?;
        let chroms = Pool::read(&mut br)?;
        let bubble_names = Pool::read(&mut br)?;
        let dbvar_ids = Pool::read(&mut br)?;

        if samples.strings.len() as u32 != n_samples
            || classes.strings.len() as u32 != n_classes
            || chroms.strings.len() as u32 != n_chroms
        {
            bail!("pool count mismatch with header");
        }

        let mut chrom_index = Vec::with_capacity(n_chroms as usize);
        for _ in 0..n_chroms {
            chrom_index.push(ChromIdx {
                chrom_name_idx: br.read_u32::<LittleEndian>()?,
                length_bp: br.read_u32::<LittleEndian>()?,
                bubble_offset: br.read_u64::<LittleEndian>()?,
                bubble_count: br.read_u32::<LittleEndian>()?,
            });
        }

        let mut bubbles = Vec::with_capacity(n_bubbles as usize);
        for _ in 0..n_bubbles {
            let start = br.read_i32::<LittleEndian>()?;
            let end = br.read_i32::<LittleEndian>()?;
            let vaf_q16 = br.read_u16::<LittleEndian>()?;
            let chrom_idx = br.read_u16::<LittleEndian>()?;
            let length_log10_q = br.read_u8()?;
            let bflags = br.read_u8()?;
            let sample_idx = br.read_u8()?;
            let class_idx = br.read_u8()?;
            let n_alts = br.read_u8()?;
            let total_reads_log = br.read_u8()?;
            let dbvar_recip_q8 = br.read_u8()?;
            let _pad = br.read_u8()?;
            let bubble_name_idx = br.read_u32::<LittleEndian>()?;
            let dbvar_id_idx = br.read_u32::<LittleEndian>()?;
            let _reserved = br.read_u32::<LittleEndian>()?;
            bubbles.push(BubbleRec {
                start, end, vaf_q16, chrom_idx,
                length_log10_q, flags: bflags, sample_idx, class_idx,
                n_alts, total_reads_log, dbvar_recip_q8,
                bubble_name_idx, dbvar_id_idx,
            });
        }

        Ok(Bbf {
            header, samples, classes, chroms,
            bubble_names, dbvar_ids,
            chrom_index, bubbles,
        })
    }

    /// Bubbles overlapping [start, end] on chrom.  Sorted by start, so we
    /// upper-bound on start <= end and scan backwards until bubble.end is
    /// also < start.  Caller still needs to filter records with end < start.
    pub fn query<'a>(&'a self, chrom_idx: u16, start: i32, end: i32) -> &'a [BubbleRec] {
        let cidx = chrom_idx as usize;
        if cidx >= self.chrom_index.len() { return &[]; }
        let ci = &self.chrom_index[cidx];
        let off = ci.bubble_offset as usize;
        let cnt = ci.bubble_count as usize;
        if cnt == 0 { return &[]; }
        let slice = &self.bubbles[off..off + cnt];

        // Upper bound: first record whose start > end.  Records >= ub can be
        // ignored entirely (their start lies after the viewport).
        let ub = slice.partition_point(|b| b.start <= end);

        // Walk back from `ub` until we find a bubble whose end is past `start`,
        // capped so a single freak wide bubble doesn't make us scan to chrom 0.
        // Since most bubbles in IGH are < 200 kb, this rarely walks more than
        // a few hundred records even at deep zoom.
        let mut lb = ub;
        let cap = ub.min(2048);
        for i in 0..cap {
            let idx = ub - 1 - i;
            if slice[idx].end < start {
                lb = idx + 1;
                break;
            }
            lb = idx;
            if i + 1 == cap { break; }
        }
        &slice[lb..ub]
    }
}
