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
    /// Population count: number of nssv observations sharing the matched
    /// REGIONID in dbVar. Proxy for "how many independent samples carry
    /// this SV". Saturates at 65535. 0 = no dbVar match (or singleton).
    pub dbvar_pop_count: u16,
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
    /// Number of dbVar nssv observations sharing this region (population
    /// frequency proxy). 0 = no dbVar match.
    pub fn pop_count(&self) -> u16 { self.dbvar_pop_count }
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
            let dbvar_pop_count = br.read_u16::<LittleEndian>()?;
            let _reserved = br.read_u16::<LittleEndian>()?;
            bubbles.push(BubbleRec {
                start, end, vaf_q16, chrom_idx,
                length_log10_q, flags: bflags, sample_idx, class_idx,
                n_alts, total_reads_log, dbvar_recip_q8,
                bubble_name_idx, dbvar_id_idx,
                dbvar_pop_count,
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn make_pool(strings: &[&str]) -> Vec<u8> {
        let mut block = Vec::new();
        let mut offsets = Vec::with_capacity(strings.len());
        for s in strings {
            offsets.push(block.len() as u32);
            block.extend_from_slice(s.as_bytes());
            block.push(0);
        }
        let mut out = Vec::new();
        out.extend_from_slice(&(block.len() as u32).to_le_bytes());
        out.extend_from_slice(&(strings.len() as u32).to_le_bytes());
        for off in &offsets {
            out.extend_from_slice(&off.to_le_bytes());
        }
        out.extend_from_slice(&block);
        out
    }

    fn make_header(n_chroms: u32, n_samples: u32, n_classes: u32, n_bubbles: u64) -> Vec<u8> {
        let mut h = Vec::with_capacity(HEADER_SZ);
        h.extend_from_slice(MAGIC);
        h.extend_from_slice(&1u16.to_le_bytes()); // version
        h.extend_from_slice(&0u16.to_le_bytes()); // flags (no zstd)
        h.extend_from_slice(&n_chroms.to_le_bytes());
        h.extend_from_slice(&n_samples.to_le_bytes());
        h.extend_from_slice(&n_classes.to_le_bytes());
        h.extend_from_slice(&n_bubbles.to_le_bytes());
        h.extend_from_slice(&0u64.to_le_bytes()); // built_unix
        let mut ref_bytes = [0u8; 16];
        ref_bytes[..6].copy_from_slice(b"GRCh38");
        h.extend_from_slice(&ref_bytes);
        h.extend_from_slice(&[0u8; 12]);
        assert_eq!(h.len(), HEADER_SZ);
        h
    }

    fn push_chrom_index(out: &mut Vec<u8>, name_idx: u32, length: u32, off: u64, count: u32) {
        out.extend_from_slice(&name_idx.to_le_bytes());
        out.extend_from_slice(&length.to_le_bytes());
        out.extend_from_slice(&off.to_le_bytes());
        out.extend_from_slice(&count.to_le_bytes());
    }

    fn push_bubble(
        out: &mut Vec<u8>,
        start: i32, end: i32,
        vaf_q16: u16, chrom_idx: u16,
        sample_idx: u8, class_idx: u8,
        bubble_name_idx: u32, dbvar_id_idx: u32, pop_count: u16,
    ) {
        out.extend_from_slice(&start.to_le_bytes());
        out.extend_from_slice(&end.to_le_bytes());
        out.extend_from_slice(&vaf_q16.to_le_bytes());
        out.extend_from_slice(&chrom_idx.to_le_bytes());
        out.push(0); // length_log10_q
        out.push(0); // flags
        out.push(sample_idx);
        out.push(class_idx);
        out.push(2); // n_alts
        out.push(0); // total_reads_log
        out.push(0); // dbvar_recip_q8
        out.push(0); // pad
        out.extend_from_slice(&bubble_name_idx.to_le_bytes());
        out.extend_from_slice(&dbvar_id_idx.to_le_bytes());
        out.extend_from_slice(&pop_count.to_le_bytes());
        out.extend_from_slice(&0u16.to_le_bytes()); // reserved
    }

    fn build_minimal_vpf() -> Vec<u8> {
        // 1 chrom, 2 samples, 3 classes, 4 bubbles
        let mut body = Vec::new();
        body.extend_from_slice(&make_pool(&["HG002", "HG00733"]));
        body.extend_from_slice(&make_pool(&["DUP", "DEL", "INS"]));
        body.extend_from_slice(&make_pool(&["chr14"]));
        body.extend_from_slice(&make_pool(&["b1", "b2", "b3", "b4"]));
        body.extend_from_slice(&make_pool(&["nssv1", ""])); // empty string at idx 1
        push_chrom_index(&mut body, 0, 107_043_718, 0, 4);
        // Bubbles must be sorted by start within a chrom.
        push_bubble(&mut body, 100,  500,  6553,  0, 0, 0, 0, 0, 0); // VAF ~0.10
        push_bubble(&mut body, 200,  800,  3276,  0, 1, 1, 1, NIL_IDX, 12);
        push_bubble(&mut body, 1000, 1500, 32767, 0, 0, 2, 2, 1, 0);  // VAF ~0.50
        push_bubble(&mut body, 2000, 2300, 65535, 0, 1, 0, 3, NIL_IDX, 200); // VAF=1.0
        let mut buf = make_header(1, 2, 3, 4);
        buf.extend_from_slice(&body);
        buf
    }

    #[test]
    fn parse_minimal_roundtrip() {
        let buf = build_minimal_vpf();
        let bbf = Bbf::parse(&buf).expect("parse");
        assert_eq!(bbf.header.n_bubbles, 4);
        assert_eq!(bbf.header.n_samples, 2);
        assert_eq!(bbf.header.n_classes, 3);
        assert_eq!(bbf.header.reference_id, "GRCh38");
        assert_eq!(bbf.samples.get(0), Some("HG002"));
        assert_eq!(bbf.samples.get(1), Some("HG00733"));
        assert_eq!(bbf.classes.get(2), Some("INS"));
        assert_eq!(bbf.chroms.get(0), Some("chr14"));
        assert_eq!(bbf.bubble_names.get(3), Some("b4"));
        assert_eq!(bbf.dbvar_ids.get(0), Some("nssv1"));
        assert_eq!(bbf.dbvar_ids.get(NIL_IDX), None);
        assert_eq!(bbf.bubbles.len(), 4);
        assert!((bbf.bubbles[2].vaf() - 0.5).abs() < 0.01);
        assert_eq!(bbf.bubbles[3].pop_count(), 200);
    }

    #[test]
    fn query_window_returns_overlapping_only() {
        let buf = build_minimal_vpf();
        let bbf = Bbf::parse(&buf).expect("parse");
        // Window [400, 1200] should overlap b1 (100-500), b2 (200-800), b3 (1000-1500); skip b4.
        let hits = bbf.query(0, 400, 1200);
        assert_eq!(hits.len(), 3);
        // Window past everything → empty.
        assert_eq!(bbf.query(0, 5000, 6000).len(), 0);
        // Out-of-range chrom_idx → empty, no panic.
        assert_eq!(bbf.query(7, 0, 100).len(), 0);
    }

    #[test]
    fn rejects_bad_magic() {
        let mut buf = build_minimal_vpf();
        buf[0] = b'X';
        assert!(Bbf::parse(&buf).is_err());
    }

    #[test]
    fn rejects_truncated_header() {
        let buf = vec![0u8; 32];
        assert!(Bbf::parse(&buf).is_err());
    }

    #[test]
    fn opens_bundled_sample() {
        let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("samples").join("igh_3sample.vpf");
        if !p.exists() {
            eprintln!("skipping: {} not present", p.display());
            return;
        }
        let bbf = Bbf::open(&p).expect("bundled sample should load");
        assert!(bbf.header.n_bubbles > 0, "sample has zero bubbles");
        assert!(bbf.header.n_samples > 0, "sample has zero samples");
        // Sanity: at least one bubble belongs to chr14.
        assert!(
            bbf.chroms.strings.iter().any(|s| s == "chr14"),
            "sample expected to contain chr14 (IGH)"
        );
    }
}
