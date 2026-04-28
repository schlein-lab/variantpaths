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
        if &magic != MAGIC { bail!("bad bbz magic: {:?}", magic); }
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
