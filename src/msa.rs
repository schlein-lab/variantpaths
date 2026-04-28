//! MSA (multi-sequence alignment) viewer state + rendering.
//!
//! When the user right-clicks "Show MSA at this locus" we collect every
//! bubble whose entry/exit anchors fall in the same 5 kb bucket as the
//! clicked bubble, pull all their .bbz ALT sequences, and pop up a window
//! with them stacked as rows.
//!
//! Two alignment modes:
//!  - PadRight (default, instant): right-pad shorter sequences with '-'.
//!    Useful as a first look; columns line up only when sequences start
//!    at the same anchor.
//!  - NeedlemanWunsch (on demand, slower): pairwise global alignment of
//!    each non-reference row against the longest row, gap-extended.
//!
//! Divergent columns are flagged in a top "*" bar.

use egui::{Color32, FontId};

#[derive(Clone)]
pub struct MsaRow {
    pub label: String,           // bubble_name + "/altN"
    pub bubble_name: String,
    pub alt_idx: u8,
    pub class: String,
    pub vaf: f32,
    /// Aligned sequence (may include '-' gaps).  Same length across rows.
    pub aligned: Vec<u8>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum AlignMode { PadRight, NeedlemanWunsch }

#[derive(Clone)]
pub struct MsaState {
    pub title: String,
    pub rows: Vec<MsaRow>,
    pub mode: AlignMode,
    /// Horizontal scroll offset (in characters) — large MSA wraps don't
    /// fit in egui's window, so we offer a slider.
    pub h_scroll_chars: usize,
    pub view_chars: usize,    // how many chars per row visible at once
}

impl MsaState {
    pub fn new<I: IntoIterator<Item = MsaRow>>(title: String, rows: I) -> Self {
        let mut rows: Vec<MsaRow> = rows.into_iter().collect();
        // initial alignment = pad-right
        let max = rows.iter().map(|r| r.aligned.len()).max().unwrap_or(0);
        for r in &mut rows {
            while r.aligned.len() < max { r.aligned.push(b'-'); }
        }
        Self {
            title,
            rows,
            mode: AlignMode::PadRight,
            h_scroll_chars: 0,
            view_chars: 200,
        }
    }

    pub fn n_cols(&self) -> usize {
        self.rows.iter().map(|r| r.aligned.len()).max().unwrap_or(0)
    }

    /// Per-column majority base + whether the column is "divergent"
    /// (i.e. not all rows agree, ignoring gaps for the consensus call).
    pub fn column_summary(&self, col: usize) -> (u8, bool) {
        let mut counts = [0u32; 256];
        let mut total: u32 = 0;
        for r in &self.rows {
            if let Some(&b) = r.aligned.get(col) {
                counts[b as usize] += 1;
                if b != b'-' { total += 1; }
            }
        }
        let (mut maj_b, mut maj_n) = (b'-', 0u32);
        for (b, n) in counts.iter().enumerate() {
            if *n > maj_n {
                maj_n = *n; maj_b = b as u8;
            }
        }
        let divergent = if total == 0 {
            false
        } else {
            // any non-majority non-gap?
            self.rows.iter().any(|r| match r.aligned.get(col) {
                Some(&b) if b != b'-' => b != maj_b,
                _ => false,
            })
        };
        (maj_b, divergent)
    }

    /// In-place re-alignment with simple Needleman-Wunsch against the
    /// longest input row, computed from the *original* (un-padded) text.
    /// Recover original by stripping '-' from current aligned form, then
    /// re-aligning each non-reference row to the reference using NW with
    /// linear gap penalty.
    pub fn realign_nw(&mut self) {
        if self.rows.len() < 2 { return; }
        // recover original sequences (strip gaps)
        let originals: Vec<Vec<u8>> = self.rows.iter()
            .map(|r| r.aligned.iter().copied().filter(|&b| b != b'-').collect())
            .collect();
        // pick the longest as reference
        let (ref_idx, _) = originals.iter().enumerate()
            .max_by_key(|(_, s)| s.len()).unwrap();

        let ref_seq = &originals[ref_idx];
        let mut new_rows = Vec::with_capacity(self.rows.len());

        // Reference row gets aligned with itself (no gaps).
        let mut ref_aligned = ref_seq.clone();

        // For NW we'll build per-row insertions of gaps in the reference;
        // then merge them into a single reference layout.
        let mut alignments: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(self.rows.len());
        for (i, seq) in originals.iter().enumerate() {
            if i == ref_idx {
                alignments.push((ref_seq.clone(), ref_seq.clone()));
            } else {
                let (ra, qa) = needleman_wunsch(ref_seq, seq);
                alignments.push((ra, qa));
            }
        }

        // Merge gap insertions across all alignments to produce a unified
        // reference layout.  This is a simplification (proper progressive
        // MSA would re-align iteratively), but it's good enough for the
        // small N (<10) we expect here.
        let merged_ref = merge_reference_gaps(&alignments);
        // Now project each query alignment onto the merged reference.
        for (i, r) in self.rows.iter().enumerate() {
            let (ra, qa) = &alignments[i];
            let projected = project_to_merged(&merged_ref, ra, qa);
            new_rows.push(MsaRow {
                label: r.label.clone(),
                bubble_name: r.bubble_name.clone(),
                alt_idx: r.alt_idx,
                class: r.class.clone(),
                vaf: r.vaf,
                aligned: projected,
            });
        }
        ref_aligned = merged_ref;
        let _ = ref_aligned;
        self.rows = new_rows;
        self.mode = AlignMode::NeedlemanWunsch;
    }
}

// ---------- NW implementation (linear gap, match=+2, mismatch=-1, gap=-2) ----------

fn needleman_wunsch(a: &[u8], b: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let n = a.len(); let m = b.len();
    if n == 0 { return (vec![b'-'; m], b.to_vec()); }
    if m == 0 { return (a.to_vec(), vec![b'-'; n]); }
    // Cap to avoid quadratic blowup on huge sequences.
    const MAX_DIM: usize = 20_000;
    if n > MAX_DIM || m > MAX_DIM {
        // Fall back to pad-right alignment.
        let max = n.max(m);
        let mut aa = a.to_vec(); aa.resize(max, b'-');
        let mut bb = b.to_vec(); bb.resize(max, b'-');
        return (aa, bb);
    }
    let cols = m + 1;
    let mut score = vec![0i32; (n + 1) * cols];
    let mut trace = vec![0u8; (n + 1) * cols];   // 0=stop 1=diag 2=up 3=left

    for i in 0..=n { score[i * cols] = -2 * i as i32; trace[i * cols] = 2; }
    for j in 0..=m { score[j] = -2 * j as i32; trace[j] = 3; }
    trace[0] = 0;

    for i in 1..=n {
        for j in 1..=m {
            let m_score = score[(i-1) * cols + (j-1)]
                + if a[i-1] == b[j-1] { 2 } else { -1 };
            let up_score = score[(i-1) * cols + j] - 2;
            let lf_score = score[i * cols + (j-1)] - 2;
            if m_score >= up_score && m_score >= lf_score {
                score[i*cols + j] = m_score; trace[i*cols + j] = 1;
            } else if up_score >= lf_score {
                score[i*cols + j] = up_score; trace[i*cols + j] = 2;
            } else {
                score[i*cols + j] = lf_score; trace[i*cols + j] = 3;
            }
        }
    }
    let mut aa = Vec::with_capacity(n + m);
    let mut bb = Vec::with_capacity(n + m);
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        match trace[i*cols + j] {
            1 => { aa.push(a[i-1]); bb.push(b[j-1]); i -= 1; j -= 1; }
            2 => { aa.push(a[i-1]); bb.push(b'-'); i -= 1; }
            3 => { aa.push(b'-'); bb.push(b[j-1]); j -= 1; }
            _ => break,
        }
    }
    aa.reverse(); bb.reverse();
    (aa, bb)
}

/// Given multiple (ref_aligned, query_aligned) pairs, merge gap insertions
/// in the reference so all share a common layout.
fn merge_reference_gaps(alignments: &[(Vec<u8>, Vec<u8>)]) -> Vec<u8> {
    // For each non-gap reference base position k, compute the maximum
    // number of gap-inserts before it across all alignments.
    // (Reference bases are assumed identical — they all derive from the
    // same original.)
    let n_ref = alignments[0].0.iter().filter(|&&b| b != b'-').count();
    let mut max_inserts = vec![0usize; n_ref + 1];
    for (ra, _) in alignments {
        let mut k = 0;
        let mut cur_inserts = 0;
        for &b in ra {
            if b == b'-' { cur_inserts += 1; }
            else {
                if cur_inserts > max_inserts[k] { max_inserts[k] = cur_inserts; }
                cur_inserts = 0;
                k += 1;
            }
        }
        if cur_inserts > max_inserts[k] { max_inserts[k] = cur_inserts; }
    }
    // Rebuild the reference layout = original ref base sequence + inserted gaps
    let ref_orig: Vec<u8> = alignments[0].0.iter().filter(|&&b| b != b'-').copied().collect();
    let mut merged = Vec::with_capacity(ref_orig.len() + max_inserts.iter().sum::<usize>());
    for k in 0..=n_ref {
        for _ in 0..max_inserts[k] { merged.push(b'-'); }
        if k < n_ref { merged.push(ref_orig[k]); }
    }
    merged
}

fn project_to_merged(merged_ref: &[u8], ra: &[u8], qa: &[u8]) -> Vec<u8> {
    // Walk merged_ref and ra in parallel, output qa chars where merged_ref
    // matches ra positions; emit '-' wherever merged_ref has a gap not in ra.
    let mut out = Vec::with_capacity(merged_ref.len());
    let mut ir = 0;
    for &m in merged_ref {
        if m == b'-' {
            // a gap in the merged ref; if ra also has a gap here we use qa, else '-'
            if ir < ra.len() && ra[ir] == b'-' {
                out.push(qa.get(ir).copied().unwrap_or(b'-'));
                ir += 1;
            } else {
                out.push(b'-');
            }
        } else {
            // walk ra's gaps until we hit a non-gap
            while ir < ra.len() && ra[ir] == b'-' {
                out.push(qa.get(ir).copied().unwrap_or(b'-'));
                ir += 1;
            }
            if ir < ra.len() {
                out.push(qa[ir]);
                ir += 1;
            }
        }
    }
    out
}

// ---------- Rendering ----------

pub fn render_msa_window(
    ctx: &egui::Context,
    state: &mut MsaState,
    open: &mut bool,
) {
    let title = state.title.clone();
    egui::Window::new(format!("MSA — {}", title))
        .open(open)
        .default_size([1100.0, 480.0])
        .resizable(true)
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label(format!("{} rows · {} cols", state.rows.len(), state.n_cols()));
                ui.separator();
                if ui.add_enabled(state.mode == AlignMode::PadRight,
                    egui::Button::new("Realign (Needleman-Wunsch)")).clicked()
                {
                    state.realign_nw();
                }
                if ui.button("Reset (pad-right)").clicked() {
                    let raw: Vec<MsaRow> = state.rows.iter().map(|r| MsaRow {
                        label: r.label.clone(),
                        bubble_name: r.bubble_name.clone(),
                        alt_idx: r.alt_idx,
                        class: r.class.clone(),
                        vaf: r.vaf,
                        aligned: r.aligned.iter().copied().filter(|&b| b != b'-').collect(),
                    }).collect();
                    *state = MsaState::new(state.title.clone(), raw);
                }
                if ui.button("Copy MSA as FASTA").clicked() {
                    let mut s = String::new();
                    for r in &state.rows {
                        s.push_str(&format!(">{}\n", r.label));
                        for chunk in r.aligned.chunks(60) {
                            s.push_str(std::str::from_utf8(chunk).unwrap_or(""));
                            s.push('\n');
                        }
                    }
                    ui.output_mut(|o| o.copied_text = s);
                }
            });
            ui.separator();

            let cw = 7.5_f32; // monospace cell width approximation
            let row_h = 14.0_f32;
            let n_cols = state.n_cols();
            let labels_w = 200.0;
            let avail_w = ui.available_width() - labels_w;
            let max_visible_cols = ((avail_w / cw) as usize).max(40).min(n_cols.max(40));

            // horizontal scroll
            ui.horizontal(|ui| {
                ui.label("scroll:");
                let mut sc = state.h_scroll_chars as i32;
                let max_sc = (n_cols as i32 - max_visible_cols as i32).max(0);
                ui.add(egui::Slider::new(&mut sc, 0..=max_sc.max(1)).integer());
                state.h_scroll_chars = (sc as usize).min(max_sc as usize);
                ui.label(format!("col {}/{}", state.h_scroll_chars, n_cols));
            });

            let start = state.h_scroll_chars;
            let end = (start + max_visible_cols).min(n_cols);

            // Divergence bar at top
            let (rect, _) = ui.allocate_exact_size(
                egui::Vec2::new(labels_w + (end - start) as f32 * cw, row_h),
                egui::Sense::hover());
            let p = ui.painter_at(rect);
            p.text(
                egui::Pos2::new(rect.left() + 4.0, rect.center().y),
                egui::Align2::LEFT_CENTER,
                "diverg.",
                FontId::monospace(10.0),
                Color32::from_gray(160));
            for col in start..end {
                let (_, div) = state.column_summary(col);
                if div {
                    let x = rect.left() + labels_w + (col - start) as f32 * cw;
                    p.line_segment(
                        [egui::Pos2::new(x + cw * 0.5, rect.top() + 1.0),
                         egui::Pos2::new(x + cw * 0.5, rect.bottom() - 1.0)],
                        egui::Stroke::new(2.0, Color32::from_rgb(0xff, 0x88, 0x44)));
                }
            }

            // Render each row
            for r in &state.rows {
                let (rect, _) = ui.allocate_exact_size(
                    egui::Vec2::new(labels_w + (end - start) as f32 * cw, row_h),
                    egui::Sense::hover());
                let p = ui.painter_at(rect);
                let lbl = format!("{}  vaf={:.2}%", r.label, r.vaf * 100.0);
                p.text(
                    egui::Pos2::new(rect.left() + 4.0, rect.center().y),
                    egui::Align2::LEFT_CENTER,
                    &lbl,
                    FontId::monospace(10.0),
                    Color32::from_gray(220));
                for col in start..end {
                    let b = r.aligned.get(col).copied().unwrap_or(b'-');
                    let (maj, _) = state.column_summary(col);
                    let bg = match b {
                        b'A' => Color32::from_rgb(0x66, 0xc2, 0xa5),
                        b'T' => Color32::from_rgb(0xfc, 0x8d, 0x62),
                        b'C' => Color32::from_rgb(0x8d, 0xa0, 0xcb),
                        b'G' => Color32::from_rgb(0xe7, 0x8a, 0xc3),
                        b'-' => Color32::from_gray(60),
                        _ =>    Color32::from_gray(110),
                    };
                    let mismatch_emphasis = b != b'-' && b != maj;
                    let x = rect.left() + labels_w + (col - start) as f32 * cw;
                    let cell = egui::Rect::from_min_size(
                        egui::Pos2::new(x, rect.top() + 1.0),
                        egui::Vec2::new(cw, row_h - 2.0));
                    let fill = if mismatch_emphasis { bg } else { bg.linear_multiply(0.35) };
                    p.rect_filled(cell, 0.0, fill);
                    if cw >= 6.0 {
                        p.text(
                            cell.center(),
                            egui::Align2::CENTER_CENTER,
                            std::str::from_utf8(&[b]).unwrap_or(" "),
                            FontId::monospace(9.0),
                            if mismatch_emphasis { Color32::WHITE }
                            else { Color32::from_gray(230) });
                    }
                }
            }
        });
}
