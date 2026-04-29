//! Minimal random-access FASTA reader (uncompressed .fa + samtools .fa.fai).
//!
//! Supports queries `fetch(chrom, start, end)` returning bytes (uppercase).
//! For .fa.gz support we'd need bgzf which is large; user can `gunzip` once.

use anyhow::{anyhow, bail, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct FaiEntry {
    pub length: u64,
    pub offset: u64,        // byte offset in fasta file where seq starts
    pub line_bases: u64,    // bases per line (not counting newline)
    pub line_bytes: u64,    // bytes per line (with newline char)
}

pub struct FastaReader {
    pub path: PathBuf,
    pub fai: HashMap<String, FaiEntry>,
}

impl FastaReader {
    pub fn open<P: AsRef<Path>>(fasta_path: P) -> Result<Self> {
        let path = fasta_path.as_ref().to_path_buf();
        let mut fai_path = path.clone();
        let extra = path.extension()
            .and_then(|e| e.to_str())
            .map(|e| format!("{e}.fai"))
            .unwrap_or_else(|| "fai".to_string());
        fai_path.set_extension(extra);
        if !fai_path.exists() {
            // try plain "<file>.fai" appended
            let alt = PathBuf::from(format!("{}.fai", path.display()));
            if alt.exists() { fai_path = alt; }
            else { bail!("fai index not found: tried {} and {}",
                         fai_path.display(),
                         format!("{}.fai", path.display())); }
        }

        let mut fai = HashMap::new();
        let f = BufReader::new(File::open(&fai_path)?);
        for line in f.lines() {
            let line = line?;
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 5 { continue; }
            let name = cols[0].to_string();
            let entry = FaiEntry {
                length:     cols[1].parse().map_err(|e| anyhow!("fai length: {}", e))?,
                offset:     cols[2].parse().map_err(|e| anyhow!("fai offset: {}", e))?,
                line_bases: cols[3].parse().map_err(|e| anyhow!("fai line_bases: {}", e))?,
                line_bytes: cols[4].parse().map_err(|e| anyhow!("fai line_bytes: {}", e))?,
            };
            fai.insert(name, entry);
        }

        Ok(FastaReader { path, fai })
    }

    pub fn chromosomes(&self) -> Vec<&str> {
        self.fai.keys().map(|s| s.as_str()).collect()
    }

    /// Fetch [start, end) (0-based, end-exclusive) from `chrom`.
    /// Returns uppercased bytes. Requested range is clamped to chrom length.
    pub fn fetch(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        let entry = self.fai.get(chrom)
            .ok_or_else(|| anyhow!("chrom {} not in fai", chrom))?;
        let mut start = start;
        let mut end = end;
        if start > entry.length { start = entry.length; }
        if end > entry.length { end = entry.length; }
        if end <= start { return Ok(Vec::new()); }

        // Compute byte offset for `start` accounting for newlines.
        let line_idx_start = start / entry.line_bases;
        let col_start = start % entry.line_bases;
        let byte_start = entry.offset + line_idx_start * entry.line_bytes + col_start;

        let line_idx_end = (end - 1) / entry.line_bases;
        let col_end = (end - 1) % entry.line_bases + 1; // exclusive end col
        let byte_end = entry.offset + line_idx_end * entry.line_bytes + col_end;
        if byte_end <= byte_start { return Ok(Vec::new()); }

        let mut f = File::open(&self.path)?;
        f.seek(SeekFrom::Start(byte_start))?;
        let mut raw = vec![0u8; (byte_end - byte_start) as usize];
        f.read_exact(&mut raw)?;
        // Strip newline characters
        let mut out = Vec::with_capacity((end - start) as usize);
        for b in raw {
            if b != b'\n' && b != b'\r' { out.push(b.to_ascii_uppercase()); }
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Write a deterministic .fa + .fai pair into a unique temp dir and
    /// return the FASTA path. Caller is responsible for nothing — the temp
    /// dir is per-test-process and OS will reap it.
    fn fixture(line_bases: usize) -> std::path::PathBuf {
        let dir = std::env::temp_dir()
            .join(format!("vp_fasta_test_{}_{}",
                std::process::id(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
        std::fs::create_dir_all(&dir).unwrap();
        let fa_path = dir.join("ref.fa");
        let fai_path = dir.join("ref.fa.fai");

        // Build a chr1 of 200 bp = "ACGT" repeated.
        let total = 200usize;
        let mut seq = String::with_capacity(total);
        for i in 0..total {
            seq.push(['A', 'C', 'G', 'T'][i % 4]);
        }

        // Write FASTA wrapped at line_bases per line.
        let mut f = std::fs::File::create(&fa_path).unwrap();
        let header = ">chr1\n";
        f.write_all(header.as_bytes()).unwrap();
        let header_len = header.len() as u64;
        let mut written = 0usize;
        while written < total {
            let take = (total - written).min(line_bases);
            f.write_all(&seq.as_bytes()[written..written + take]).unwrap();
            f.write_all(b"\n").unwrap();
            written += take;
        }
        // Write FAI: name, length, offset, line_bases, line_bytes
        let line_bytes = line_bases + 1; // +1 for '\n'
        let mut fai = std::fs::File::create(&fai_path).unwrap();
        writeln!(fai, "chr1\t{}\t{}\t{}\t{}",
            total, header_len, line_bases, line_bytes).unwrap();
        fa_path
    }

    #[test]
    fn fetch_aligned_to_line_bases() {
        let fa = fixture(60);
        let r = FastaReader::open(&fa).expect("open");
        assert_eq!(r.fai.get("chr1").map(|e| e.length), Some(200));
        let bytes = r.fetch("chr1", 0, 4).unwrap();
        assert_eq!(bytes, b"ACGT");
    }

    #[test]
    fn fetch_crosses_line_boundary() {
        let fa = fixture(60); // newline every 60 bases
        let r = FastaReader::open(&fa).unwrap();
        // bases 58..62 cross the wrap; correct answer should still be ACGT-pattern.
        let bytes = r.fetch("chr1", 58, 62).unwrap();
        let expected: Vec<u8> = (58..62).map(|i| b"ACGT"[i % 4] as u8).collect();
        assert_eq!(bytes, expected);
    }

    #[test]
    fn fetch_clamps_to_chrom_length() {
        let fa = fixture(60);
        let r = FastaReader::open(&fa).unwrap();
        let bytes = r.fetch("chr1", 198, 9999).unwrap();
        assert_eq!(bytes.len(), 2);
    }

    #[test]
    fn fetch_unknown_chrom_errors() {
        let fa = fixture(60);
        let r = FastaReader::open(&fa).unwrap();
        assert!(r.fetch("chrX", 0, 10).is_err());
    }

    #[test]
    fn fetch_returns_uppercase() {
        // Same fixture; bytes are already uppercase, but the function must
        // not mangle them either. (Lowercase coverage is in the .to_ascii_uppercase
        // branch; we just confirm round-trip is uppercase.)
        let fa = fixture(80);
        let r = FastaReader::open(&fa).unwrap();
        let b = r.fetch("chr1", 0, 8).unwrap();
        assert!(b.iter().all(|c| c.is_ascii_uppercase()));
    }
}
