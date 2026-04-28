//! Application-level model: loaded data + selection state.

use std::collections::HashMap;

use crate::format::{annot::AnnotTrack, bbf::Bbf, bbz::Bbz, fasta::FastaReader};

const DENSITY_REF_BIN_BP: i32 = 5_000;

/// All loaded data.  Anything optional may be `None`.
pub struct Atlas {
    pub bbf: Option<Bbf>,
    pub bbf_path: Option<String>,
    pub bbz: Option<Bbz>,
    pub bbz_path: Option<String>,
    pub reference: Option<FastaReader>,
    pub reference_path: Option<String>,
    pub annot_tracks: Vec<AnnotTrack>,

    /// Stable y-axis upper bound per (chrom_idx, sample_idx) pair, in
    /// "bubbles per `DENSITY_REF_BIN_BP` (5 kb) bin".  Computed once at
    /// load.  Density-mode rendering scales bar heights against this so
    /// the y-axis does not auto-rescale when the user pans away from a
    /// dense region.
    pub density_max_per_5kb: HashMap<(u16, u8), u32>,
}

impl Atlas {
    pub fn empty() -> Self {
        Self {
            bbf: None,
            bbf_path: None,
            bbz: None,
            bbz_path: None,
            reference: None,
            reference_path: None,
            annot_tracks: vec![crate::format::annot::igh_default_track()],
            density_max_per_5kb: HashMap::new(),
        }
    }

    /// Recompute the per-(chrom, sample) max-bin density.  Call after
    /// loading or replacing the BBF.
    pub fn recompute_density_scale(&mut self) {
        self.density_max_per_5kb.clear();
        let Some(bbf) = self.bbf.as_ref() else { return };
        // bin counters keyed by (chrom_idx, sample_idx, bin)
        let mut counts: HashMap<(u16, u8, i32), u32> = HashMap::new();
        for b in &bbf.bubbles {
            let mid = (b.start as i64 + b.end as i64) / 2;
            let bin = (mid / DENSITY_REF_BIN_BP as i64) as i32;
            *counts.entry((b.chrom_idx, b.sample_idx, bin)).or_insert(0) += 1;
        }
        for ((c, s, _), n) in counts {
            let e = self.density_max_per_5kb.entry((c, s)).or_insert(0);
            if n > *e { *e = n; }
        }
    }

    /// Returns the y-axis upper bound at the *current* bin size, scaled
    /// from the stored 5kb-resolution max.  When the renderer's bin size
    /// differs from 5 kb we extrapolate proportionally so peaks stay at
    /// roughly full track height across all zoom levels.
    pub fn density_y_max(&self, chrom_idx: u16, sample_idx: u8, bin_bp: f64) -> f32 {
        let m = *self.density_max_per_5kb.get(&(chrom_idx, sample_idx)).unwrap_or(&1);
        let scale = (bin_bp / DENSITY_REF_BIN_BP as f64).max(1.0) as f32;
        ((m as f32) * scale).max(1.0)
    }

    /// 5kb-resolution max (for label display).
    pub fn density_ref_max(&self, chrom_idx: u16, sample_idx: u8) -> u32 {
        *self.density_max_per_5kb.get(&(chrom_idx, sample_idx)).unwrap_or(&1)
    }

    /// (analysis_start_bp, analysis_end_bp) for a chrom: the leftmost
    /// bubble.start and rightmost bubble.end in the BBF.  Anything
    /// outside this range is "unmapped" — no BRANCH analysis happened
    /// there.  Distinct from `length_bp` in the chrom_index, which is
    /// just the max(end).
    pub fn analysis_range(&self, chrom_idx: u16) -> Option<(i32, i32)> {
        let bbf = self.bbf.as_ref()?;
        let ci = bbf.chrom_index.get(chrom_idx as usize)?;
        let off = ci.bubble_offset as usize;
        let cnt = ci.bubble_count as usize;
        if cnt == 0 { return None; }
        let slice = &bbf.bubbles[off..off + cnt];
        let start_min = slice.iter().map(|b| b.start).min()?;
        let end_max = slice.iter().map(|b| b.end).max()?;
        Some((start_min, end_max))
    }

    /// Map a class *name* to a stable color (independent of pool order).
    pub fn class_color_by_name(&self, name: &str) -> egui::Color32 {
        match name {
            "GERMLINE_KNOWN"          => egui::Color32::from_rgb(0x2c, 0x7f, 0xb8), // dark blue
            "GERMLINE_NOVEL"          => egui::Color32::from_rgb(0x7f, 0xcd, 0xbb), // teal
            "SUBCLONAL_KNOWN"         => egui::Color32::from_rgb(0xfd, 0xae, 0x61), // orange
            "SUBCLONAL_NOVEL"         => egui::Color32::from_rgb(0x98, 0x4e, 0xa3), // purple
            "UNCLASSIFIED"            => egui::Color32::from_rgb(0x99, 0x99, 0x99),
            _ => {
                // deterministic fallback hash for unknown classes
                let mut h: u32 = 2166136261;
                for &b in name.as_bytes() { h = h.wrapping_mul(16777619) ^ b as u32; }
                let r = ( h        & 0xFF) as u8;
                let g = ((h >>  8) & 0xFF) as u8;
                let b = ((h >> 16) & 0xFF) as u8;
                egui::Color32::from_rgb(r, g, b)
            }
        }
    }

    /// Look up the class name in the BBF pool and color it.
    pub fn class_color(&self, class_idx: u8) -> egui::Color32 {
        let name = self.bbf.as_ref()
            .and_then(|b| b.classes.strings.get(class_idx as usize))
            .cloned()
            .unwrap_or_default();
        self.class_color_by_name(&name)
    }
}
