//! All canvas rendering. One function per layer.
//!
//! Coordinate system:
//!   bp_to_x(bp) = canvas.left + (bp - view_start) / span * canvas.width
//!
//! The View's bp_per_px controls which bubble layer is active.

use egui::{Align2, Color32, FontId, Pos2, Rect, Stroke, Vec2};

use crate::format::bbf::{Bbf, BubbleRec};
use crate::format::fasta::FastaReader;
use crate::format::annot::AnnotTrack;
use crate::model::Atlas;
use crate::view::ViewState;

pub const RULER_H: f32 = 22.0;
pub const SEQ_H: f32 = 16.0;
pub const ANNOT_H: f32 = 22.0;
pub const SAMPLE_H: f32 = 180.0;
pub const SAMPLE_LABEL_W: f32 = 130.0;
pub const LANE_H: f32 = 12.0;
pub const LANE_GAP: f32 = 2.0;

/// bp -> screen x  (returns None if span is 0).
#[inline]
pub fn bp_to_x(bp: f64, view: &ViewState, track_rect: Rect) -> f32 {
    let frac = (bp - view.view_start_bp) / view.span();
    track_rect.left() + (frac as f32) * track_rect.width()
}
#[inline]
pub fn x_to_bp(x: f32, view: &ViewState, track_rect: Rect) -> f64 {
    let frac = ((x - track_rect.left()) / track_rect.width().max(1.0)) as f64;
    view.view_start_bp + frac * view.span()
}

// ------------------- Unmapped overlay -------------------

/// Mark genome regions outside the analyzed range.  Pan keeps working
/// there, but a dark dim and a label make clear that "no bubbles" means
/// "analysis ended", not "no data".
pub fn draw_unmapped_overlay(p: &egui::Painter, rect: Rect, label: &str) {
    if rect.width() < 1.0 { return; }
    // Translucent dark fill — tracks remain dimly visible underneath.
    p.rect_filled(rect, 0.0,
        Color32::from_rgba_premultiplied(30, 30, 30, 180));
    // Single soft border at the analysis-edge side so the boundary is
    // unambiguous without diagonal noise across the bubble area.
    if rect.width() > 60.0 {
        p.text(
            rect.center(),
            Align2::CENTER_CENTER,
            label,
            FontId::proportional(11.0),
            Color32::from_rgb(220, 200, 100),
        );
    }
}

// ------------------- Ruler -------------------

pub fn draw_ruler(ui: &egui::Painter, rect: Rect, view: &ViewState) {
    ui.rect_filled(rect, 0.0, Color32::from_gray(28));
    let span = view.span();
    // Pick a tick step that gives ~6-12 ticks across the view.
    let target_ticks = 8.0_f64;
    let raw_step = span / target_ticks;
    let pow = 10f64.powf(raw_step.log10().floor());
    let step = if raw_step / pow < 2.0 { 1.0 * pow }
               else if raw_step / pow < 5.0 { 2.0 * pow }
               else { 5.0 * pow };
    let first = (view.view_start_bp / step).ceil() * step;
    let mut t = first;
    while t < view.view_end_bp {
        let x = bp_to_x(t, view, rect);
        ui.line_segment(
            [Pos2::new(x, rect.bottom() - 6.0), Pos2::new(x, rect.bottom())],
            Stroke::new(1.0, Color32::from_gray(160)),
        );
        let label = format_bp(t as i64);
        ui.text(
            Pos2::new(x + 2.0, rect.top() + 2.0),
            Align2::LEFT_TOP,
            label,
            FontId::monospace(10.0),
            Color32::from_gray(200),
        );
        t += step;
    }
}

pub fn format_bp(bp: i64) -> String {
    if bp.abs() >= 1_000_000 {
        format!("{:.2} Mb", (bp as f64) / 1e6)
    } else if bp.abs() >= 1_000 {
        format!("{:.1} kb", (bp as f64) / 1e3)
    } else {
        format!("{} bp", bp)
    }
}

// ------------------- Sequence -------------------

pub fn draw_sequence(
    ui: &egui::Painter, rect: Rect, view: &ViewState,
    fa: &FastaReader, chrom_name: &str,
) {
    ui.rect_filled(rect, 0.0, Color32::from_gray(20));
    let bp_per_px = view.bp_per_px(rect.width());
    if bp_per_px > 5.0 {
        // too zoomed out
        ui.text(
            rect.center(),
            Align2::CENTER_CENTER,
            format!("(zoom in below 5 bp/px to see sequence — current {:.1})", bp_per_px),
            FontId::proportional(10.0),
            Color32::from_gray(120),
        );
        return;
    }
    let start = view.view_start_bp.floor().max(0.0) as u64;
    let end = view.view_end_bp.ceil().max(start as f64 + 1.0) as u64;
    let seq = match fa.fetch(chrom_name, start, end) {
        Ok(s) => s,
        Err(e) => {
            ui.text(
                rect.center(),
                Align2::CENTER_CENTER,
                format!("[seq error: {}]", e),
                FontId::proportional(10.0),
                Color32::from_rgb(200, 80, 80),
            );
            return;
        }
    };
    let px_per_bp = 1.0 / (bp_per_px as f32);
    for (i, b) in seq.iter().enumerate() {
        let bp = start as f64 + i as f64;
        let x = bp_to_x(bp, view, rect);
        let col = match *b {
            b'A' => Color32::from_rgb(0x66, 0xc2, 0xa5),
            b'T' => Color32::from_rgb(0xfc, 0x8d, 0x62),
            b'C' => Color32::from_rgb(0x8d, 0xa0, 0xcb),
            b'G' => Color32::from_rgb(0xe7, 0x8a, 0xc3),
            _ =>    Color32::from_gray(110),
        };
        if px_per_bp >= 6.0 {
            // wide enough to draw a letter
            ui.text(
                Pos2::new(x + px_per_bp / 2.0, rect.center().y),
                Align2::CENTER_CENTER,
                std::str::from_utf8(&[*b]).unwrap_or("?"),
                FontId::monospace(10.0),
                col,
            );
        } else {
            // just a colored tick
            ui.line_segment(
                [Pos2::new(x, rect.top() + 2.0), Pos2::new(x, rect.bottom() - 2.0)],
                Stroke::new(px_per_bp.max(0.5), col),
            );
        }
    }
}

// ------------------- Annotation track -------------------

pub fn draw_annotation(ui: &egui::Painter, rect: Rect, view: &ViewState,
    chrom_name: &str, track: &AnnotTrack)
{
    ui.rect_filled(rect, 0.0, Color32::from_gray(24));
    let features = track.query(chrom_name, view.view_start_bp as i32, view.view_end_bp as i32);
    let cy = rect.center().y;
    for f in features {
        let x0 = bp_to_x(f.start as f64, view, rect).max(rect.left());
        let x1 = bp_to_x(f.end as f64, view, rect).min(rect.right());
        if x1 - x0 < 1.0 { continue; }
        let h = (rect.height() * 0.5).max(6.0);
        let bar = Rect::from_min_size(Pos2::new(x0, cy - h/2.0), Vec2::new(x1 - x0, h));
        ui.rect_filled(bar, 2.0, Color32::from_rgb(0x4a, 0x90, 0xe2));
        // arrow strand marker
        let arrow = if f.strand == '-' { "◀" } else { "▶" };
        let label = format!("{} {}", f.name, arrow);
        if x1 - x0 > 40.0 {
            ui.text(
                Pos2::new(x0 + 4.0, cy),
                Align2::LEFT_CENTER,
                label,
                FontId::proportional(10.0),
                Color32::WHITE,
            );
        }
    }
}

// ------------------- Bubble Track (auto-mode) -------------------

pub struct BubbleHit {
    pub global_index: usize,
    pub distance_px: f32,
}

pub fn draw_bubble_track(
    p: &egui::Painter, rect: Rect, view: &ViewState, atlas: &Atlas,
    bbf: &Bbf, sample_idx: u8, hover_pos: Option<Pos2>, mode: BubbleMode,
    selected: Option<usize>,
) -> Option<BubbleHit> {
    p.rect_filled(rect, 0.0, Color32::from_gray(18));
    p.line_segment([Pos2::new(rect.left(), rect.bottom()),
                    Pos2::new(rect.right(), rect.bottom())],
                   Stroke::new(0.5, Color32::from_gray(60)));

    // No more VAF y-axis grid: in stack mode the vertical dimension
    // encodes lane (overlap), not VAF.  VAF is shown as text on each box.

    // Slice of bubbles intersecting view.
    let bubbles = bbf.query(
        view.chrom_idx,
        view.view_start_bp as i32,
        view.view_end_bp as i32,
    );
    let _chrom_offset = bbf.chrom_index.get(view.chrom_idx as usize)
        .map(|c| c.bubble_offset as usize)
        .unwrap_or(0);

    // we need original index per bubble; we'll rebuild it as we iterate
    // (the slice is contiguous starting at chrom_offset minus walk_back, so
    // the displayed first bubble's global idx == chrom_offset + (bubble_pos
    // in the chrom slice for that record).  Rather than track the offset
    // shift, we just compute the global index from the BubbleRec address.)
    let all_bubbles_ptr = bbf.bubbles.as_ptr();
    let mut hit: Option<BubbleHit> = None;

    match mode {
        BubbleMode::Density => {
            // Bin into ~rect.width() / 6 columns.
            let n_bins = ((rect.width() / 6.0) as usize).max(8).min(2000);
            let bin_bp = view.span() / n_bins as f64;
            let mut bin_total = vec![0u32; n_bins];
            let n_classes = atlas.bbf.as_ref().map(|b| b.classes.strings.len()).unwrap_or(6);
            let mut bin_class = vec![0u32; n_bins * n_classes];
            for b in bubbles {
                if b.sample_idx != sample_idx { continue; }
                if !passes(view, b) { continue; }
                let mid = (b.start as f64 + b.end as f64) * 0.5;
                let frac = (mid - view.view_start_bp) / view.span();
                if !(0.0..=1.0).contains(&frac) { continue; }
                let bi = ((frac * n_bins as f64).floor() as usize).min(n_bins - 1);
                bin_total[bi] += 1;
                bin_class[bi * n_classes + b.class_idx as usize] += 1;
            }
            // Stable y-axis: chromosome-wide max (per sample) at 5kb resolution,
            // scaled to the current bin size.  This means the y-axis does NOT
            // change as the user pans — small peaks no longer get squashed by
            // distant high peaks like the IGHM cluster.
            let y_max = atlas.density_y_max(view.chrom_idx, sample_idx, bin_bp);
            let bin_w = rect.width() / n_bins as f32;
            let track_h = rect.height() - 4.0;
            for (bi, &tot) in bin_total.iter().enumerate() {
                if tot == 0 { continue; }
                let x0 = rect.left() + bi as f32 * bin_w;
                let mut y = rect.bottom();
                // sqrt scaling: small peaks remain visible alongside large ones
                let height_total = ((tot as f32 / y_max).clamp(0.0, 1.0)).sqrt() * track_h;
                for ci in 0..n_classes {
                    let c = bin_class[bi * n_classes + ci];
                    if c == 0 { continue; }
                    let h = (c as f32 / tot as f32) * height_total;
                    let bar = Rect::from_min_size(
                        Pos2::new(x0, y - h),
                        Vec2::new(bin_w.max(1.0) - 0.5, h));
                    p.rect_filled(bar, 0.0, atlas.class_color(ci as u8));
                    y -= h;
                }
            }
            // Y-axis label
            let max_5kb = atlas.density_ref_max(view.chrom_idx, sample_idx);
            let label = format!("n / {:.0} kb (max {} per 5kb, sqrt-scaled)",
                                bin_bp / 1000.0, max_5kb);
            p.text(
                Pos2::new(rect.right() - 4.0, rect.top() + 2.0),
                Align2::RIGHT_TOP,
                label,
                FontId::monospace(9.0),
                Color32::from_gray(160),
            );
        }
        BubbleMode::Stack | BubbleMode::Detail => {
            // ---- Haplotype baseline ----
            // The reference haplotype as a continuous horizontal line at
            // the bottom of the track.  Branches/bubbles depart from it
            // upwards into stacked lanes.
            let baseline_y = rect.bottom() - 2.0;
            p.line_segment(
                [Pos2::new(rect.left(), baseline_y),
                 Pos2::new(rect.right(), baseline_y)],
                Stroke::new(1.5, Color32::from_rgb(0x55, 0xaa, 0xff)));
            p.text(
                Pos2::new(rect.right() - 6.0, baseline_y - 1.0),
                Align2::RIGHT_BOTTOM,
                "ref haplotype",
                FontId::monospace(8.0),
                Color32::from_rgb(0x55, 0xaa, 0xff));

            // ---- Greedy lane assignment, left-to-right ----
            let n_lanes_max = ((rect.height() - 8.0) / (LANE_H + LANE_GAP)) as usize;
            let mut lane_last_x: Vec<f32> = Vec::with_capacity(n_lanes_max);
            let mut overflow: u32 = 0;
            let detail_mode = matches!(mode, BubbleMode::Detail);

            for b in bubbles {
                if b.sample_idx != sample_idx { continue; }
                if !passes(view, b) { continue; }
                let x0 = bp_to_x(b.start as f64, view, rect);
                let x1 = bp_to_x(b.end as f64, view, rect);
                // Cull only sub-pixel-AND-off-screen.  Anything that has
                // any horizontal footprint we want to draw.
                if x1 < rect.left() - 4.0 || x0 > rect.right() + 4.0 { continue; }
                let xw = (x1 - x0).max(2.0);

                // Find earliest free lane (binary lane occupancy).
                let lane = lane_last_x.iter().position(|&lx| lx + 2.0 < x0);
                let lane_idx = if let Some(li) = lane {
                    lane_last_x[li] = x1;
                    li
                } else if lane_last_x.len() < n_lanes_max {
                    lane_last_x.push(x1);
                    lane_last_x.len() - 1
                } else {
                    overflow += 1;
                    continue;
                };

                let y_top = baseline_y - 6.0 - (lane_idx + 1) as f32 * (LANE_H + LANE_GAP);
                let y_bot = y_top + LANE_H;
                let col = atlas.class_color(b.class_idx);

                // Class-name lookup (shared with germline detection)
                let cls_name = bbf.classes.get(b.class_idx as u32).unwrap_or("");
                let is_germline = cls_name.starts_with("GERMLINE");

                if is_germline {
                    // Two parallel lines = haplotype-pair representation
                    let y_a = y_top + 2.0;
                    let y_b = y_bot - 2.0;
                    p.line_segment([Pos2::new(x0, y_a), Pos2::new(x1, y_a)],
                                   Stroke::new(1.6, col));
                    p.line_segment([Pos2::new(x0, y_b), Pos2::new(x1, y_b)],
                                   Stroke::new(1.6, col));
                } else {
                    // Subclonal/unclassified: filled box
                    let box_rect = Rect::from_min_size(
                        Pos2::new(x0, y_top), Vec2::new(xw, LANE_H));
                    p.rect_filled(box_rect, 1.5, col.linear_multiply(0.55));
                    p.rect_stroke(box_rect, 1.5, Stroke::new(1.0, col));
                }

                // Branch connector: from baseline up to this lane.
                let mid = (x0 + x1) * 0.5;
                p.line_segment(
                    [Pos2::new(mid, baseline_y),
                     Pos2::new(mid, y_bot)],
                    Stroke::new(0.7, col.linear_multiply(0.45)));

                // Anchor ticks at exact start/end on baseline
                p.line_segment(
                    [Pos2::new(x0, baseline_y - 3.0),
                     Pos2::new(x0, baseline_y + 1.0)],
                    Stroke::new(1.0, col));
                p.line_segment(
                    [Pos2::new(x1, baseline_y - 3.0),
                     Pos2::new(x1, baseline_y + 1.0)],
                    Stroke::new(1.0, col));

                // Shared marker (a small dot at the right edge of the bar)
                if b.is_shared() && xw > 6.0 {
                    p.circle_filled(Pos2::new(x1 - 3.0, (y_top + y_bot) * 0.5),
                                    2.0,
                                    Color32::from_rgb(0xff, 0xff, 0x80));
                }

                // VAF label inside the box if there's room.
                if xw > 32.0 {
                    let recip_lbl = if b.dbvar_recip() > 0.0 {
                        format!("  dbVar-match {:.0}%", b.dbvar_recip() * 100.0)
                    } else { String::new() };
                    // Make it explicit that this is reciprocal-topology overlap,
                    // NOT population allele frequency.
                    let label = format!("VAF {:.2}%{}", b.vaf() * 100.0, recip_lbl);
                    p.text(
                        Pos2::new(x0 + 3.0, (y_top + y_bot) * 0.5),
                        Align2::LEFT_CENTER,
                        label,
                        FontId::monospace(9.0),
                        Color32::WHITE,
                    );
                }
                if detail_mode && xw > 80.0 {
                    if let Some(name) = bbf.bubble_names.get(b.bubble_name_idx) {
                        p.text(
                            Pos2::new(x1 - 3.0, (y_top + y_bot) * 0.5),
                            Align2::RIGHT_CENTER,
                            name,
                            FontId::monospace(8.0),
                            Color32::from_gray(220),
                        );
                    }
                }

                // Selection highlight
                if is_global_selected(b, all_bubbles_ptr, selected) {
                    let r = Rect::from_min_size(
                        Pos2::new(x0 - 1.0, y_top - 1.0),
                        Vec2::new(xw + 2.0, LANE_H + 2.0));
                    p.rect_stroke(r, 2.5, Stroke::new(2.0, Color32::WHITE));
                }

                // Hit-test (rectangular)
                if let Some(hp) = hover_pos {
                    let r = Rect::from_min_size(
                        Pos2::new(x0, y_top), Vec2::new(xw, LANE_H));
                    if r.expand(2.0).contains(hp) {
                        let gidx = ((b as *const BubbleRec) as usize
                                  - all_bubbles_ptr as usize)
                                  / std::mem::size_of::<BubbleRec>();
                        hit = Some(BubbleHit { global_index: gidx, distance_px: 0.0 });
                    }
                }
            }

            if overflow > 0 {
                p.text(
                    Pos2::new(rect.right() - 4.0, rect.top() + 2.0),
                    Align2::RIGHT_TOP,
                    format!("+{} overflow (track full)", overflow),
                    FontId::monospace(9.0),
                    Color32::from_rgb(0xff, 0x80, 0x80),
                );
            }
        }
    }

    let _ = selected;
    hit
}

fn is_global_selected(b: &BubbleRec, all_ptr: *const BubbleRec, sel: Option<usize>) -> bool {
    match sel {
        Some(idx) => {
            let g = ((b as *const BubbleRec) as usize - all_ptr as usize)
                  / std::mem::size_of::<BubbleRec>();
            g == idx
        }
        None => false,
    }
}

/// Apply ViewState filters to a single bubble.  Used in all three modes.
#[inline]
fn passes(view: &ViewState, b: &BubbleRec) -> bool {
    view.passes_filters(b)
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum BubbleMode { Density, Stack, Detail }

pub fn pick_mode(bp_per_px: f64) -> BubbleMode {
    // Density when truly zoomed out (whole chromosome / genome).
    // Stack as the working mode (each bubble in its own lane).
    // Detail is identical to Stack with bp-level text + sequence hints.
    if bp_per_px > 2000.0 { BubbleMode::Density }
    else if bp_per_px > 0.5 { BubbleMode::Stack }
    else { BubbleMode::Detail }
}

fn draw_arc(p: &egui::Painter, x0: f32, x1: f32, baseline_y: f32, height: f32,
            color: Color32, selected: bool)
{
    // Same 24-segment polyline as before — visually identical to the line-
    // segment version.  Difference: we emit it as a *single* Shape::line
    // (one tessellation pass, one PathShape in egui's mesh) instead of 24
    // independent line_segment calls.  No change in pixels rendered.
    let n = 24usize;
    let mut points = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let t = i as f32 / n as f32;
        let x = x0 + (x1 - x0) * t;
        let lift = -4.0 * height * t * (1.0 - t);
        points.push(Pos2::new(x, baseline_y + lift));
    }
    let stroke = if selected {
        Stroke::new(2.5, Color32::WHITE)
    } else {
        Stroke::new(1.4, color)
    };
    p.add(egui::Shape::line(points, stroke));
}
