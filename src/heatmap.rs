//! Cross-sample bubble-density heatmap.
//!
//! A compact strip rendered above the per-sample lane tracks: one row per
//! sample, x = genomic position (matched to the current view), color = log-
//! scaled count of bubbles passing the user's filters in that bin.  This
//! makes "which samples carry more variation in this region" answerable at
//! a glance, without scrolling through the per-sample stacks.
//!
//! Notes:
//! - Color ramp: dark navy → blue → magenta → yellow (viridis-ish 4-stop).
//!   Distinct from the per-sample class colors so the heatmap reads as a
//!   different layer, not "more bubbles".
//! - Log scaling so a 1000× hot region doesn't drown a 10× warm region.
//! - Bins are computed per-frame from the same `bbf.query` window the rest
//!   of the renderer uses, so the heatmap stays in sync with pan/zoom.

use egui::{Align2, Color32, FontId, Pos2, Rect, Stroke, Vec2};

use crate::format::bbf::Bbf;
use crate::view::ViewState;

pub const SAMPLE_LABEL_W: f32 = 110.0;
pub const ROW_H: f32 = 14.0;
pub const ROW_GAP: f32 = 1.0;
pub const HEADER_H: f32 = 18.0;

/// Total pixel height needed to render `n_samples` rows.  Caller uses this
/// when allocating space in the central panel.
pub fn height(n_samples: usize) -> f32 {
    HEADER_H + n_samples as f32 * (ROW_H + ROW_GAP) + 4.0
}

/// Render the heatmap into `rect`.  No-ops for an empty atlas.
pub fn draw_heatmap(p: &egui::Painter, rect: Rect, view: &ViewState, bbf: &Bbf) {
    p.rect_filled(rect, 0.0, Color32::from_gray(14));
    p.line_segment(
        [Pos2::new(rect.left(), rect.bottom()), Pos2::new(rect.right(), rect.bottom())],
        Stroke::new(0.5, Color32::from_gray(60)),
    );

    let n_samples = bbf.header.n_samples as usize;
    if n_samples == 0 || rect.width() < SAMPLE_LABEL_W + 20.0 { return; }

    let track_left = rect.left() + SAMPLE_LABEL_W;
    let track_w = rect.width() - SAMPLE_LABEL_W;
    let n_bins = ((track_w / 4.0) as usize).clamp(8, 2000);

    // Per-(sample, bin) counts.  Keep it as a flat vec for cache.
    let mut bins = vec![0u32; n_samples * n_bins];
    let mut max_count: u32 = 1;

    let span = view.span().max(1.0);
    let bubbles = bbf.query(
        view.chrom_idx,
        view.view_start_bp as i32,
        view.view_end_bp as i32,
    );
    for b in bubbles {
        if !view.passes_filters(b) { continue; }
        let s = b.sample_idx as usize;
        if s >= n_samples { continue; }
        let mid = (b.start as f64 + b.end as f64) * 0.5;
        let frac = (mid - view.view_start_bp) / span;
        if !(0.0..=1.0).contains(&frac) { continue; }
        let bi = ((frac * n_bins as f64).floor() as usize).min(n_bins - 1);
        let idx = s * n_bins + bi;
        bins[idx] = bins[idx].saturating_add(1);
        if bins[idx] > max_count { max_count = bins[idx]; }
    }

    // ----- Header row: title + max count + bin width -----
    let title = format!(
        "cross-sample density — max {}/bin · {:.0} bp/bin",
        max_count,
        span / n_bins as f64,
    );
    p.text(
        Pos2::new(rect.left() + 6.0, rect.top() + HEADER_H * 0.5),
        Align2::LEFT_CENTER,
        title,
        FontId::monospace(10.0),
        Color32::from_gray(180),
    );
    // Tiny color-ramp legend at the right edge of the header.
    let legend_w = 80.0;
    let legend_x = rect.right() - legend_w - 6.0;
    let legend_y = rect.top() + HEADER_H * 0.5 - 4.0;
    for i in 0..(legend_w as usize) {
        let t = i as f32 / legend_w;
        let col = ramp(t);
        p.line_segment(
            [Pos2::new(legend_x + i as f32, legend_y),
             Pos2::new(legend_x + i as f32, legend_y + 8.0)],
            Stroke::new(1.0, col),
        );
    }
    p.text(Pos2::new(legend_x - 4.0, legend_y + 4.0),
        Align2::RIGHT_CENTER, "0",
        FontId::monospace(9.0), Color32::from_gray(140));
    p.text(Pos2::new(legend_x + legend_w + 4.0, legend_y + 4.0),
        Align2::LEFT_CENTER, format!("{}", max_count),
        FontId::monospace(9.0), Color32::from_gray(140));

    // ----- Sample rows -----
    let bin_w = track_w / n_bins as f32;
    let log_max = ((max_count as f32).ln() + 1.0).max(1.0);

    for s in 0..n_samples {
        let y0 = rect.top() + HEADER_H + s as f32 * (ROW_H + ROW_GAP);
        let y1 = y0 + ROW_H;
        // Sample label
        let name = bbf.samples.get(s as u32).unwrap_or("?");
        p.text(
            Pos2::new(rect.left() + 6.0, (y0 + y1) * 0.5),
            Align2::LEFT_CENTER,
            name,
            FontId::monospace(10.0),
            Color32::from_gray(210),
        );
        // Track background
        let track_rect = Rect::from_min_size(
            Pos2::new(track_left, y0), Vec2::new(track_w, y1 - y0));
        p.rect_filled(track_rect, 0.0, Color32::from_gray(22));
        // Bins
        for bi in 0..n_bins {
            let c = bins[s * n_bins + bi];
            if c == 0 { continue; }
            let intensity = ((c as f32).ln() + 1.0) / log_max;
            let col = ramp(intensity.clamp(0.0, 1.0));
            let x0 = track_left + bi as f32 * bin_w;
            let bar = Rect::from_min_size(
                Pos2::new(x0, y0), Vec2::new(bin_w.max(1.0), ROW_H));
            p.rect_filled(bar, 0.0, col);
        }
    }
}

/// 4-stop perceptual-ish ramp.  Not viridis exactly — picked so the dark
/// end blends with the panel background, the bright end is unambiguous,
/// and the midtones don't compete with the per-class category colors.
fn ramp(t: f32) -> Color32 {
    const STOPS: &[(f32, [u8; 3])] = &[
        (0.00, [12, 16, 40]),
        (0.40, [40, 60, 180]),
        (0.70, [200, 60, 160]),
        (1.00, [240, 220, 80]),
    ];
    let t = t.clamp(0.0, 1.0);
    // Find segment.
    let mut i = 0;
    while i + 1 < STOPS.len() && t > STOPS[i + 1].0 { i += 1; }
    let (t0, c0) = STOPS[i];
    let (t1, c1) = STOPS[(i + 1).min(STOPS.len() - 1)];
    let span = (t1 - t0).max(1e-6);
    let f = ((t - t0) / span).clamp(0.0, 1.0);
    let lerp = |a: u8, b: u8| (a as f32 * (1.0 - f) + b as f32 * f).round() as u8;
    Color32::from_rgb(
        lerp(c0[0], c1[0]),
        lerp(c0[1], c1[1]),
        lerp(c0[2], c1[2]),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ramp_endpoints_are_stable() {
        let a = ramp(0.0);
        let b = ramp(1.0);
        // Dark end should be very dark.
        assert!(a.r() < 40 && a.g() < 40 && a.b() < 80, "{a:?}");
        // Hot end should be bright yellow-ish.
        assert!(b.r() > 200 && b.g() > 180 && b.b() < 120, "{b:?}");
    }

    #[test]
    fn ramp_clamps_out_of_range() {
        let lo = ramp(-1.0);
        let hi = ramp(2.0);
        assert_eq!((lo.r(), lo.g(), lo.b()), (ramp(0.0).r(), ramp(0.0).g(), ramp(0.0).b()));
        assert_eq!((hi.r(), hi.g(), hi.b()), (ramp(1.0).r(), ramp(1.0).g(), ramp(1.0).b()));
    }

    #[test]
    fn height_scales_with_samples() {
        assert!(height(0) < height(3));
        assert!(height(3) < height(10));
        // A reasonable cap so we don't try to draw a 10k-pixel strip.
        assert!(height(50) < 1000.0);
    }
}
