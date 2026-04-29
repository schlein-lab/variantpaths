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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_bed(contents: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!("vp_annot_{}_{}.bed",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(contents.as_bytes()).unwrap();
        path
    }

    #[test]
    fn parses_minimal_bed() {
        let p = write_bed("\
chr1\t100\t200\tfeatA\t.\t+
chr1\t300\t400\tfeatB\t.\t-
# comment line
chr2\t50\t150\tfeatC\t.\t+
");
        let t = AnnotTrack::open_bed(&p).expect("open");
        assert_eq!(t.features.len(), 3);
        // Sorted by chrom then start.
        assert_eq!(t.features[0].name, "featA");
        assert_eq!(t.features[1].name, "featB");
        assert_eq!(t.features[2].chrom, "chr2");
    }

    #[test]
    fn skips_track_and_browser_headers() {
        let p = write_bed("\
track name=foo
browser hide all
chr3\t10\t20\tx\t.\t+
");
        let t = AnnotTrack::open_bed(&p).unwrap();
        assert_eq!(t.features.len(), 1);
    }

    #[test]
    fn three_column_bed_uses_default_name_and_strand() {
        let p = write_bed("chr1\t10\t20\n");
        let t = AnnotTrack::open_bed(&p).unwrap();
        assert_eq!(t.features.len(), 1);
        assert_eq!(t.features[0].name, ".");
        assert_eq!(t.features[0].strand, '.');
    }

    #[test]
    fn query_returns_overlapping_only() {
        let p = write_bed("\
chr1\t100\t200\tA
chr1\t300\t400\tB
chr1\t500\t600\tC
");
        let t = AnnotTrack::open_bed(&p).unwrap();
        let hits = t.query("chr1", 150, 350);
        let names: Vec<&str> = hits.iter().map(|f| f.name.as_str()).collect();
        assert_eq!(names, vec!["A", "B"]);
        assert_eq!(t.query("chr2", 0, 1_000_000).len(), 0);
        // Half-open: feature ending exactly at start is excluded.
        assert_eq!(t.query("chr1", 200, 250).len(), 0);
    }

    #[test]
    fn igh_default_track_has_expected_genes() {
        let t = igh_default_track();
        let names: Vec<&str> = t.features.iter().map(|f| f.name.as_str()).collect();
        for g in &["IGHG4", "IGHG2", "IGHG1", "IGHG3", "IGHM", "IGHD"] {
            assert!(names.contains(g), "missing {} in built-in track", g);
        }
        assert!(t.features.iter().all(|f| f.chrom == "chr14"));
    }
}
