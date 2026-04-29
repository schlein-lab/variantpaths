//! .vpz reader (VariantPaths sequence Z-archive). Companion to .bbf.
//!
//! See build_vpz.py for the writer.  Format:
//!   Header (32 byte LE):
//!     magic            "VPZ1"
//!     version          u16 = 1
//!     flags            u16  bit0 = zstd payload
//!     n_bubbles        u32
//!     bbf_built_unix   u64
//!     reserved         u32 * 3
//!   Body (zstd-decompressed):
//!     for each bubble:
//!       u32 name_len; bytes name
//!       u8  n_alts
//!       for each alt:
//!         u32 path_len_nodes; u32[path_len] node_ids
//!         u32 alt_seq_len; bytes alt_seq

use anyhow::{anyhow, bail, Result};
use byteorder::{LittleEndian, ReadBytesExt};
use std::collections::HashMap;
use std::fs::File;
use std::io::{Cursor, Read};
use std::path::Path;

pub const MAGIC: &[u8; 4] = b"VPZ1";
/// Legacy magic from the prototype phase, before the on-disk extension was
/// renamed from .bbz to .vpz. Files written by the pre-rename build script
/// still carry these four bytes; we accept them as long as the body layout
/// is identical, which it has been since version 1.
pub const MAGIC_LEGACY: &[u8; 4] = b"BBZ1";
pub const HEADER_SZ: usize = 32;
pub const FLAG_ZSTD: u16 = 1 << 0;

#[derive(Debug, Clone)]
pub struct AltSeq {
    pub path_nodes: Vec<u32>,
    pub seq: Vec<u8>,
}

#[derive(Debug, Clone)]
pub struct Bbz {
    pub version: u16,
    pub n_bubbles: u32,
    pub bbf_built_unix: u64,
    /// Per-bubble alt sequences, keyed by bubble_name.
    pub by_name: HashMap<String, Vec<AltSeq>>,
}

impl Bbz {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut buf = Vec::new();
        File::open(path.as_ref())?.read_to_end(&mut buf)?;
        Self::parse(&buf)
    }

    pub fn parse(buf: &[u8]) -> Result<Self> {
        if buf.len() < HEADER_SZ { bail!("bbz too small for header"); }
        let mut rdr = Cursor::new(&buf[..HEADER_SZ]);
        let mut magic = [0u8; 4];
        rdr.read_exact(&mut magic)?;
        if &magic != MAGIC && &magic != MAGIC_LEGACY {
            bail!("bad bbz magic: {:?}", magic);
        }
        let version = rdr.read_u16::<LittleEndian>()?;
        let flags   = rdr.read_u16::<LittleEndian>()?;
        let n_bubbles = rdr.read_u32::<LittleEndian>()?;
        let bbf_built_unix = rdr.read_u64::<LittleEndian>()?;
        let _reserved0 = rdr.read_u32::<LittleEndian>()?;
        let _reserved1 = rdr.read_u32::<LittleEndian>()?;
        let _reserved2 = rdr.read_u32::<LittleEndian>()?;

        let body = if flags & FLAG_ZSTD != 0 {
            zstd::decode_all(&buf[HEADER_SZ..])
                .map_err(|e| anyhow!("bbz zstd decode: {}", e))?
        } else {
            buf[HEADER_SZ..].to_vec()
        };

        let mut br = Cursor::new(body.as_slice());
        let mut by_name: HashMap<String, Vec<AltSeq>> = HashMap::new();
        for _ in 0..n_bubbles {
            let name_len = br.read_u32::<LittleEndian>()? as usize;
            let mut name_buf = vec![0u8; name_len];
            br.read_exact(&mut name_buf)?;
            let name = String::from_utf8(name_buf)
                .map_err(|e| anyhow!("bbz bubble_name utf8: {}", e))?;
            let n_alts = br.read_u8()? as usize;
            let mut alts = Vec::with_capacity(n_alts);
            for _ in 0..n_alts {
                let plen = br.read_u32::<LittleEndian>()? as usize;
                let mut path_nodes = Vec::with_capacity(plen);
                for _ in 0..plen {
                    path_nodes.push(br.read_u32::<LittleEndian>()?);
                }
                let slen = br.read_u32::<LittleEndian>()? as usize;
                let mut seq = vec![0u8; slen];
                br.read_exact(&mut seq)?;
                alts.push(AltSeq { path_nodes, seq });
            }
            by_name.insert(name, alts);
        }

        Ok(Bbz {
            version, n_bubbles, bbf_built_unix, by_name,
        })
    }

    pub fn alts_for(&self, bubble_name: &str) -> Option<&[AltSeq]> {
        self.by_name.get(bubble_name).map(|v| v.as_slice())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn make_header(n_bubbles: u32) -> Vec<u8> {
        let mut h = Vec::with_capacity(HEADER_SZ);
        h.extend_from_slice(MAGIC);
        h.extend_from_slice(&1u16.to_le_bytes()); // version
        h.extend_from_slice(&0u16.to_le_bytes()); // flags (no zstd)
        h.extend_from_slice(&n_bubbles.to_le_bytes());
        h.extend_from_slice(&0u64.to_le_bytes()); // bbf_built_unix
        h.extend_from_slice(&[0u8; 12]); // reserved
        assert_eq!(h.len(), HEADER_SZ);
        h
    }

    fn push_alt(out: &mut Vec<u8>, path_nodes: &[u32], seq: &[u8]) {
        out.extend_from_slice(&(path_nodes.len() as u32).to_le_bytes());
        for n in path_nodes { out.extend_from_slice(&n.to_le_bytes()); }
        out.extend_from_slice(&(seq.len() as u32).to_le_bytes());
        out.extend_from_slice(seq);
    }

    fn push_bubble(out: &mut Vec<u8>, name: &str, alts: &[(&[u32], &[u8])]) {
        out.extend_from_slice(&(name.len() as u32).to_le_bytes());
        out.extend_from_slice(name.as_bytes());
        out.push(alts.len() as u8);
        for (path, seq) in alts { push_alt(out, path, seq); }
    }

    fn build_minimal_vpz() -> Vec<u8> {
        let mut body = Vec::new();
        push_bubble(&mut body, "b1", &[
            (&[10u32, 11, 12], b"ATCG"),
            (&[10, 13, 12],   b"ATGG"),
        ]);
        push_bubble(&mut body, "b_long_name", &[
            (&[100], b"AAAA"),
        ]);
        let mut buf = make_header(2);
        buf.extend_from_slice(&body);
        buf
    }

    #[test]
    fn parse_minimal_roundtrip() {
        let buf = build_minimal_vpz();
        let bbz = Bbz::parse(&buf).expect("parse");
        assert_eq!(bbz.n_bubbles, 2);
        let alts = bbz.alts_for("b1").expect("b1 present");
        assert_eq!(alts.len(), 2);
        assert_eq!(alts[0].path_nodes, vec![10, 11, 12]);
        assert_eq!(alts[0].seq, b"ATCG");
        assert_eq!(alts[1].seq, b"ATGG");
        let alts2 = bbz.alts_for("b_long_name").expect("b_long_name present");
        assert_eq!(alts2.len(), 1);
        assert_eq!(alts2[0].seq, b"AAAA");
        assert!(bbz.alts_for("missing").is_none());
    }

    #[test]
    fn rejects_bad_magic() {
        let mut buf = build_minimal_vpz();
        buf[0] = b'X';
        assert!(Bbz::parse(&buf).is_err());
    }

    #[test]
    fn opens_bundled_sample() {
        let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("samples").join("igh_3sample.vpz");
        if !p.exists() {
            eprintln!("skipping: {} not present", p.display());
            return;
        }
        let bbz = Bbz::open(&p).expect("bundled sample should load");
        assert!(bbz.n_bubbles > 0, "sample has zero bubbles");
        assert!(!bbz.by_name.is_empty(), "sample by_name map is empty");
    }
}
