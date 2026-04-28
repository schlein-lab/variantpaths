//! Minimal BED parser.  GFF3 deferred (BED is enough for the MVP).

use anyhow::{anyhow, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct AnnotFeature {
    pub chrom: String,
    pub start: i32,
    pub end: i32,
    pub name: String,
    pub strand: char,
}

#[derive(Default, Clone)]
pub struct AnnotTrack {
    pub label: String,
    pub features: Vec<AnnotFeature>,
}

impl AnnotTrack {
    pub fn open_bed<P: AsRef<Path>>(path: P) -> Result<Self> {
        let label = path.as_ref().file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("annot")
            .to_string();
        let f = BufReader::new(File::open(path.as_ref())?);
        let mut features = Vec::new();
        for line in f.lines() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') || line.starts_with("track")
                || line.starts_with("browser") {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 3 { continue; }
            let chrom = cols[0].to_string();
            let start: i32 = cols[1].parse().map_err(|e| anyhow!("bed start: {}", e))?;
            let end:   i32 = cols[2].parse().map_err(|e| anyhow!("bed end: {}", e))?;
            let name = cols.get(3).unwrap_or(&".").to_string();
            let strand = cols.get(5)
                .and_then(|s| s.chars().next())
                .unwrap_or('.');
            features.push(AnnotFeature { chrom, start, end, name, strand });
        }
        // Sort per chrom, by start
        features.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start)));
        Ok(AnnotTrack { label, features })
    }

    /// Features overlapping [start, end] on chrom.
    pub fn query<'a>(&'a self, chrom: &str, start: i32, end: i32) -> Vec<&'a AnnotFeature> {
        // linear scan; for a few thousand features this is fine.  If a user
        // loads a 1M-record BED we'd swap in a per-chrom sorted index.
        self.features.iter()
            .filter(|f| f.chrom == chrom && f.end > start && f.start < end)
            .collect()
    }
}

/// Hard-coded fallback IGH gene annotation.  Used when no BED is loaded.
pub fn igh_default_track() -> AnnotTrack {
    let chrom = "chr14".to_string();
    AnnotTrack {
        label: "IGH (built-in)".to_string(),
        features: vec![
            mk(&chrom, 105_586_000, 105_590_000, "IGHE",  '+'),
            mk(&chrom, 105_618_000, 105_622_000, "IGHG4", '+'),
            mk(&chrom, 105_642_000, 105_647_000, "IGHG2", '+'),
            mk(&chrom, 105_652_000, 105_656_000, "IGHA2", '+'),
            mk(&chrom, 105_711_000, 105_714_000, "IGHA1", '+'),
            mk(&chrom, 105_745_000, 105_750_000, "IGHG1", '+'),
            mk(&chrom, 105_770_000, 105_774_000, "IGHG3", '+'),
            mk(&chrom, 105_854_000, 105_864_000, "IGHM",  '+'),
            mk(&chrom, 105_867_000, 105_877_000, "IGHD",  '+'),
            mk(&chrom, 105_877_000, 106_300_000, "V/D/J", '+'),
        ],
    }
}

fn mk(chrom: &str, start: i32, end: i32, name: &str, strand: char) -> AnnotFeature {
    AnnotFeature { chrom: chrom.to_string(), start, end, name: name.to_string(), strand }
}
