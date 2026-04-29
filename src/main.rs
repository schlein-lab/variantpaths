#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod format;
mod model;
mod view;
mod render;
mod msa;

use eframe::egui;
use std::path::PathBuf;

use crate::format::{annot::AnnotTrack, bbf::Bbf, bbz::Bbz, fasta::FastaReader};
use crate::model::Atlas;
use crate::view::{SharedFilter, ViewState};
use crate::render::{
    bp_to_x, x_to_bp, draw_ruler, draw_sequence, draw_annotation,
    draw_bubble_track, draw_unmapped_overlay, pick_mode, format_bp,
    RULER_H, SEQ_H, ANNOT_H, SAMPLE_LABEL_W, LANE_H, LANE_GAP,
};

struct App {
    atlas: Atlas,
    view: ViewState,
    error_msg: Option<String>,
    info_msg: Option<String>,
    drag_active: bool,
    /// Rolling timing buffers for the on-screen perf overlay (visible
    /// when --perf or F11).
    perf_visible: bool,
    last_frame_ms: f32,
    last_query_ms: f32,
    last_render_ms: f32,
    last_visible_n: u32,
    msa_state: Option<msa::MsaState>,
    msa_open: bool,
}

impl App {
    fn new(_cc: &eframe::CreationContext<'_>, args: Vec<String>) -> Self {
        let perf = args.iter().any(|a| a == "--perf");
        let mut app = Self {
            atlas: Atlas::empty(),
            view: ViewState::fit_chrom(1),
            error_msg: None,
            info_msg: None,
            drag_active: false,
            perf_visible: perf,
            last_frame_ms: 0.0,
            last_query_ms: 0.0,
            last_render_ms: 0.0,
            last_visible_n: 0,
            msa_state: None,
            msa_open: false,
        };
        // Process CLI args by extension
        for arg in args.iter().skip(1) {
            app.open_path(PathBuf::from(arg));
        }
        app
    }

    fn open_path(&mut self, path: PathBuf) {
        let lower = path.to_string_lossy().to_lowercase();
        if lower.ends_with(".vpf") {
            self.load_bbf(&path);
        } else if lower.ends_with(".vpz") {
            self.load_bbz(&path);
        } else if lower.ends_with(".fa") || lower.ends_with(".fasta")
                  || lower.ends_with(".fna") {
            self.load_fasta(&path);
        } else if lower.ends_with(".bed") {
            self.load_bed(&path);
        } else {
            self.error_msg = Some(format!("unknown file type: {}", path.display()));
        }
    }

    fn load_bbz(&mut self, path: &PathBuf) {
        match Bbz::open(path) {
            Ok(z) => {
                self.info_msg = Some(format!(
                    "loaded {} (alt-sequences for {} bubbles)",
                    path.display(), z.n_bubbles));
                self.atlas.bbz = Some(z);
                self.atlas.bbz_path = Some(path.display().to_string());
            }
            Err(e) => self.error_msg = Some(format!("VPZ load failed: {}", e)),
        }
    }

    fn load_bbf(&mut self, path: &PathBuf) {
        match Bbf::open(path) {
            Ok(b) => {
                self.info_msg = Some(format!(
                    "loaded {} ({} bubbles, {} samples, {} chroms, ref={})",
                    path.display(),
                    b.header.n_bubbles, b.header.n_samples, b.header.n_chroms,
                    b.header.reference_id,
                ));
                let chrom_len = b.chrom_index.first().map(|c| c.length_bp).unwrap_or(1);
                // Reset all view + filter state to defaults so the previous file's
                // VAF/length/recurrence/pop-count sliders don't silently hide records
                // in the new dataset.
                self.view = ViewState::fit_chrom(chrom_len);
                self.view.sample_visible = vec![true; b.header.n_samples as usize];
                self.view.class_visible = vec![true; b.header.n_classes as usize];
                // For IGH demo: zoom to interesting range if it looks like chr14:IGH
                if b.chroms.strings.first().map(|s| s.as_str()) == Some("chr14")
                    && chrom_len > 105_000_000
                {
                    self.view.view_start_bp = 105_400_000.0;
                    self.view.view_end_bp = 106_500_000.0;
                }
                self.atlas.bbf = Some(b);
                self.atlas.bbf_path = Some(path.display().to_string());
                self.atlas.recompute_density_scale();
            }
            Err(e) => self.error_msg = Some(format!("VPF load failed: {}", e)),
        }
    }

    fn load_fasta(&mut self, path: &PathBuf) {
        match FastaReader::open(path) {
            Ok(fa) => {
                self.info_msg = Some(format!(
                    "loaded reference {} ({} chroms)",
                    path.display(), fa.fai.len()
                ));
                self.atlas.reference = Some(fa);
                self.atlas.reference_path = Some(path.display().to_string());
            }
            Err(e) => self.error_msg = Some(format!("FASTA load failed: {}", e)),
        }
    }

    fn load_bed(&mut self, path: &PathBuf) {
        match AnnotTrack::open_bed(path) {
            Ok(t) => {
                self.info_msg = Some(format!("loaded annot {} ({} features)",
                    path.display(), t.features.len()));
                self.atlas.annot_tracks.push(t);
            }
            Err(e) => self.error_msg = Some(format!("BED load failed: {}", e)),
        }
    }

    fn current_chrom_name(&self) -> Option<&str> {
        self.atlas.bbf.as_ref()
            .and_then(|b| b.chroms.strings.get(self.view.chrom_idx as usize))
            .map(|s| s.as_str())
    }

    fn current_chrom_length(&self) -> Option<u32> {
        self.atlas.bbf.as_ref()
            .and_then(|b| b.chrom_index.get(self.view.chrom_idx as usize))
            .map(|c| c.length_bp)
    }

    fn jump_to(&mut self, q: &str) {
        // accepts:  "chr14:105_678_492"  "chr14"  "105678492"
        let q = q.trim().replace('_', "").replace(',', "");
        if let Some((chrom, pos)) = q.split_once(':') {
            self.set_chrom_by_name(chrom);
            if let Ok(p) = pos.parse::<f64>() {
                self.center_at(p);
            }
        } else if let Ok(p) = q.parse::<f64>() {
            self.center_at(p);
        } else if !q.is_empty() {
            self.set_chrom_by_name(&q);
        }
    }

    fn set_chrom_by_name(&mut self, name: &str) {
        if let Some(b) = &self.atlas.bbf {
            if let Some(idx) = b.chroms.strings.iter().position(|c| c == name) {
                self.view.chrom_idx = idx as u16;
                let len = b.chrom_index.get(idx).map(|c| c.length_bp).unwrap_or(1);
                self.view.view_start_bp = 0.0;
                self.view.view_end_bp = len as f64;
            }
        }
    }

    fn center_at(&mut self, bp: f64) {
        let span = self.view.span();
        self.view.view_start_bp = bp - span / 2.0;
        self.view.view_end_bp = self.view.view_start_bp + span;
        self.view.clamp(self.current_chrom_length().map(|l| l as f64));
    }

    fn n_samples(&self) -> usize {
        self.atlas.bbf.as_ref().map(|b| b.header.n_samples as usize).unwrap_or(0)
    }

    /// Collect MSA rows for all bubbles at the same locus as the given
    /// bubble.  "Same locus" = entry/exit fall in the same 5 kb bucket.
    /// Pulls each bubble's alt-sequences from the loaded .bbz; rows that
    /// have no .bbz entry are silently skipped.
    fn open_msa_at_locus(&mut self, anchor_idx: usize) {
        let bbf = match &self.atlas.bbf { Some(b) => b, None => return };
        let bbz = match &self.atlas.bbz { Some(z) => z, None => {
            self.error_msg = Some(
                "MSA needs a .vpz loaded (Open .vpz first)".to_string()); return;
        }};
        let anchor = match bbf.bubbles.get(anchor_idx) { Some(b) => *b, None => return };
        const BUCKET: i32 = 5_000;
        let a_bin = anchor.start / BUCKET;
        let b_bin = anchor.end / BUCKET;
        let chrom = bbf.chroms.get(anchor.chrom_idx as u32).unwrap_or("?").to_string();

        let mut rows: Vec<msa::MsaRow> = Vec::new();
        for b in &bbf.bubbles {
            if b.chrom_idx != anchor.chrom_idx { continue; }
            if b.start / BUCKET != a_bin || b.end / BUCKET != b_bin { continue; }
            let bname = bbf.bubble_names.get(b.bubble_name_idx)
                .map(|s| s.to_string());
            let bname = match bname { Some(n) if !n.is_empty() => n, _ => continue };
            let alts = match bbz.alts_for(&bname) { Some(a) => a, None => continue };
            let class = bbf.classes.get(b.class_idx as u32).unwrap_or("?").to_string();
            let sample = bbf.samples.get(b.sample_idx as u32).unwrap_or("?").to_string();
            for (i, alt) in alts.iter().enumerate() {
                rows.push(msa::MsaRow {
                    label: format!("{}/{} alt{}  {}", sample, &bname, i, &class),
                    bubble_name: bname.clone(),
                    alt_idx: i as u8,
                    class: class.clone(),
                    vaf: b.vaf(),
                    aligned: alt.seq.clone(),
                });
            }
        }

        if rows.is_empty() {
            self.error_msg = Some(format!(
                "no .vpz alt-sequences at this locus ({}:{}..{})",
                chrom, anchor.start, anchor.end));
            return;
        }
        let title = format!(
            "{}:{}..{}  ({} alt-rows from {} bubbles)",
            chrom, anchor.start, anchor.end,
            rows.len(),
            rows.iter().map(|r| &r.bubble_name).collect::<std::collections::HashSet<_>>().len()
        );
        self.msa_state = Some(msa::MsaState::new(title, rows));
        self.msa_open = true;
    }
}

impl eframe::App for App {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let frame_t0 = std::time::Instant::now();

        // F11 toggles the perf overlay
        ctx.input(|i| {
            if i.key_pressed(egui::Key::F11) { self.perf_visible = !self.perf_visible; }
        });
        // ----- TopPanel: toolbar -----
        egui::TopBottomPanel::top("topbar").show(ctx, |ui| {
            // first row: file / nav
            ui.horizontal(|ui| {
                if ui.button("Open .vpf").clicked() {
                    if let Some(p) = rfd::FileDialog::new()
                        .add_filter("VariantPaths File", &["vpf"])
                        .pick_file() { self.load_bbf(&p); }
                }
                if ui.button("Open Reference (.fa)").clicked() {
                    if let Some(p) = rfd::FileDialog::new()
                        .add_filter("FASTA", &["fa","fasta","fna"])
                        .pick_file() { self.load_fasta(&p); }
                }
                if ui.button("Open BED").clicked() {
                    if let Some(p) = rfd::FileDialog::new()
                        .add_filter("BED", &["bed"])
                        .pick_file() { self.load_bed(&p); }
                }
                if ui.button("Open .vpz (sequences)").clicked() {
                    if let Some(p) = rfd::FileDialog::new()
                        .add_filter("VariantPaths Sequences", &["vpz"])
                        .pick_file() { self.load_bbz(&p); }
                }
                ui.separator();
                let has_bbf = self.atlas.bbf.is_some();
                if has_bbf {
                    let (chrom_names, chrom_lens, current_idx) = {
                        let b = self.atlas.bbf.as_ref().unwrap();
                        (
                            b.chroms.strings.clone(),
                            b.chrom_index.iter().map(|c| c.length_bp).collect::<Vec<_>>(),
                            self.view.chrom_idx,
                        )
                    };
                    let mut new_chrom_idx: Option<u16> = None;
                    egui::ComboBox::from_label("chrom")
                        .selected_text(chrom_names.get(current_idx as usize)
                            .cloned().unwrap_or_default())
                        .show_ui(ui, |ui| {
                            for (i, name) in chrom_names.iter().enumerate() {
                                if ui.selectable_label(i as u16 == current_idx, name).clicked() {
                                    new_chrom_idx = Some(i as u16);
                                }
                            }
                        });
                    if let Some(i) = new_chrom_idx {
                        self.view.chrom_idx = i;
                        let len = chrom_lens.get(i as usize).copied().unwrap_or(1);
                        self.view.view_start_bp = 0.0;
                        self.view.view_end_bp = len as f64;
                    }
                    ui.add_sized([240.0, 20.0],
                        egui::TextEdit::singleline(&mut self.view.jump_input)
                            .hint_text("e.g. chr14:105_700_000"));
                    if ui.button("Go").clicked() {
                        let s = self.view.jump_input.clone();
                        self.jump_to(&s);
                    }
                    ui.label(format!(
                        "{}:{} - {}",
                        chrom_names.get(self.view.chrom_idx as usize)
                            .map(|s| s.as_str()).unwrap_or(""),
                        format_bp(self.view.view_start_bp as i64),
                        format_bp(self.view.view_end_bp as i64),
                    ));
                } else {
                    ui.label("no .vpf loaded — drop a file or use Open .vpf");
                }
            });

            // ----- Second row: Filters -----
            if let Some(b) = self.atlas.bbf.as_ref() {
                let sample_names = b.samples.strings.clone();
                let class_names = b.classes.strings.clone();
                drop(b); // release immutable borrow

                ui.horizontal_wrapped(|ui| {
                    egui::CollapsingHeader::new("Filters")
                        .default_open(true)
                        .show(ui, |ui| {
                        ui.horizontal_wrapped(|ui| {
                            // VAF range
                            ui.label("VAF:");
                            ui.add(egui::Slider::new(&mut self.view.vaf_min, 0.0..=1.0)
                                .text("min")
                                .logarithmic(true)
                                .min_decimals(4)
                                .step_by(0.0001));
                            ui.add(egui::Slider::new(&mut self.view.vaf_max, 0.0..=1.0)
                                .text("max")
                                .logarithmic(true)
                                .min_decimals(4)
                                .step_by(0.0001));
                            if self.view.vaf_min > self.view.vaf_max {
                                self.view.vaf_max = self.view.vaf_min;
                            }
                            // VAF presets
                            if ui.small_button("subclonal <5%").clicked() {
                                self.view.vaf_min = 0.0; self.view.vaf_max = 0.05;
                            }
                            if ui.small_button("intermediate 5-25%").clicked() {
                                self.view.vaf_min = 0.05; self.view.vaf_max = 0.25;
                            }
                            if ui.small_button("germline-like ≥25%").clicked() {
                                self.view.vaf_min = 0.25; self.view.vaf_max = 1.0;
                            }
                            if ui.small_button("all").clicked() {
                                self.view.vaf_min = 0.0; self.view.vaf_max = 1.0;
                            }
                        });
                        ui.horizontal(|ui| {
                            // min total_reads (back-of-envelope: log2(R)*8)
                            ui.label("min total_reads:");
                            let mut approx_reads = if self.view.min_total_reads_log == 0 { 0u32 }
                                else { (2.0_f32).powf(self.view.min_total_reads_log as f32 / 8.0) as u32 };
                            if ui.add(egui::Slider::new(&mut approx_reads, 0..=500)
                                .text("reads"))
                                .changed()
                            {
                                self.view.min_total_reads_log = if approx_reads == 0 { 0 }
                                    else { ((approx_reads as f32).log2() * 8.0).round().clamp(0.0, 255.0) as u8 };
                            }
                            ui.separator();
                            ui.label("length range (bp):");
                            ui.add(egui::DragValue::new(&mut self.view.min_length_bp)
                                .speed(100).clamp_range(0..=1_000_000_000));
                            ui.label("–");
                            ui.add(egui::DragValue::new(&mut self.view.max_length_bp)
                                .speed(100).clamp_range(0..=i32::MAX));
                            if ui.small_button("reset").clicked() {
                                self.view.min_length_bp = 0;
                                self.view.max_length_bp = i32::MAX;
                                self.view.min_total_reads_log = 0;
                            }
                        });

                        // ensure visibility vecs are sized to current sample/class count
                        if self.view.sample_visible.len() != sample_names.len() {
                            self.view.sample_visible = vec![true; sample_names.len()];
                        }
                        if self.view.class_visible.len() != class_names.len() {
                            self.view.class_visible = vec![true; class_names.len()];
                        }

                        ui.horizontal_wrapped(|ui| {
                            ui.label("recurrence:");
                            ui.radio_value(&mut self.view.shared_filter,
                                SharedFilter::All, "all");
                            ui.radio_value(&mut self.view.shared_filter,
                                SharedFilter::SharedOnly, "shared (≥2 samples)");
                            ui.radio_value(&mut self.view.shared_filter,
                                SharedFilter::PrivateOnly, "private (1 sample)");
                        });
                        ui.horizontal_wrapped(|ui| {
                            ui.label("max lanes per sample:");
                            let mut v = self.view.max_lanes_per_sample;
                            ui.add(egui::Slider::new(&mut v, 0u32..=50)
                                .integer()
                                .text("(0 = all)"));
                            self.view.max_lanes_per_sample = v;
                            if ui.small_button("6").clicked() { self.view.max_lanes_per_sample = 6; }
                            if ui.small_button("12").clicked() { self.view.max_lanes_per_sample = 12; }
                            if ui.small_button("all").clicked() { self.view.max_lanes_per_sample = 0; }
                        });
                        ui.horizontal_wrapped(|ui| {
                            ui.label("max dbVar pop count:");
                            // u16::MAX maps to "no cap"; show as ∞.
                            let mut cap_int = if self.view.max_pop_count == u16::MAX {
                                10000u32
                            } else {
                                self.view.max_pop_count as u32
                            };
                            let resp = ui.add(egui::Slider::new(&mut cap_int, 0u32..=10000)
                                .integer()
                                .logarithmic(true)
                                .text("nssv"));
                            if resp.changed() {
                                self.view.max_pop_count = if cap_int >= 10000 { u16::MAX } else { cap_int as u16 };
                            }
                            if ui.small_button("rare ≤5").clicked() { self.view.max_pop_count = 5; }
                            if ui.small_button("≤50").clicked() { self.view.max_pop_count = 50; }
                            if ui.small_button("∞").clicked() { self.view.max_pop_count = u16::MAX; }
                        });
                        ui.horizontal_wrapped(|ui| {
                            ui.checkbox(&mut self.view.force_haplotype_mode,
                                "haplotype mode (two parallel lines, bubbles routed by VAF≷0.5)");
                        });
                        ui.horizontal_wrapped(|ui| {
                            ui.label("samples:");
                            for (i, n) in sample_names.iter().enumerate() {
                                ui.checkbox(&mut self.view.sample_visible[i], n);
                            }
                        });
                        ui.horizontal_wrapped(|ui| {
                            ui.label("classes:");
                            for (i, n) in class_names.iter().enumerate() {
                                let col = self.atlas.class_color(i as u8);
                                let label = egui::RichText::new(n).color(col);
                                ui.checkbox(&mut self.view.class_visible[i], label);
                            }
                        });
                    });
                });
            }
        });

        // ----- BottomPanel: status -----
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if let Some(idx) = self.view.hover_bubble {
                    if let Some(b) = &self.atlas.bbf {
                        if let Some(rec) = b.bubbles.get(idx) {
                            let name = b.bubble_names.get(rec.bubble_name_idx).unwrap_or("");
                            let cls = b.classes.get(rec.class_idx as u32).unwrap_or("?");
                            let smp = b.samples.get(rec.sample_idx as u32).unwrap_or("?");
                            let dbv = b.dbvar_ids.get(rec.dbvar_id_idx).unwrap_or("");
                            let pop_str = if rec.pop_count() == 0 {
                                "novel".to_string()
                            } else {
                                format!("{} nssv", rec.pop_count())
                            };
                            ui.label(format!(
                                "{} | {} | {} | start={} end={} length={} VAF={:.4} | dbVar={} recip={:.2} pop={}",
                                smp, name, cls,
                                format_bp(rec.start as i64),
                                format_bp(rec.end as i64),
                                format_bp(rec.length() as i64),
                                rec.vaf(), dbv, rec.dbvar_recip(), pop_str,
                            ));
                        }
                    }
                } else if let Some(m) = &self.info_msg {
                    ui.colored_label(egui::Color32::LIGHT_GREEN, m);
                } else if let Some(e) = &self.error_msg {
                    ui.colored_label(egui::Color32::LIGHT_RED, e);
                } else {
                    ui.label("scroll = zoom (centered on cursor) | drag = pan | hover bubble for details");
                }
            });
        });

        // ----- CentralPanel: canvas -----
        // Step 1: collect dropped files (mutable self ops)
        let dropped: Vec<PathBuf> = ctx.input(|i| {
            i.raw.dropped_files.iter()
                .filter_map(|d| d.path.clone())
                .collect()
        });
        for p in dropped { self.open_path(p); }

        // Step 2: capture interaction inputs into plain values (no borrow holds)
        let scroll_y = ctx.input(|i| i.raw_scroll_delta.y);
        let key_plus = ctx.input(|i| i.key_pressed(egui::Key::Equals) || i.key_pressed(egui::Key::Plus));
        let key_minus = ctx.input(|i| i.key_pressed(egui::Key::Minus));
        let key_left = ctx.input(|i| i.key_pressed(egui::Key::ArrowLeft));
        let key_right = ctx.input(|i| i.key_pressed(egui::Key::ArrowRight));
        let key_f = ctx.input(|i| i.key_pressed(egui::Key::F));

        egui::CentralPanel::default().show(ctx, |ui| {
            // Empty placeholder if no data
            if self.atlas.bbf.is_none() {
                ui.centered_and_justified(|ui| {
                    ui.label("No data loaded. Drop a .vpf file here or use Open .vpf.");
                });
                return;
            }

            // ---- Pre-pass: BP-space lane assignment per sample ----
            // We need this to compute per-sample track heights BEFORE we
            // allocate the canvas, so the ScrollArea can size correctly
            // and show every bubble.  Lane assignment in BP-space is
            // conservative (slightly more lanes than the pixel-space
            // assignment in the actual renderer) — that's fine because
            // it ensures we never run out of vertical room.
            let n_samples = self.atlas.bbf.as_ref().unwrap().header.n_samples as usize;
            let mut visible_per_sample: Vec<u32> = vec![0; n_samples];
            let mut total_per_sample: Vec<u32> = vec![0; n_samples];
            let mut lanes_used_per_sample: Vec<usize> = vec![0; n_samples];
            let mut total_visible_for_perf: u32 = 0;
            let query_t0 = std::time::Instant::now();
            {
                let bbf = self.atlas.bbf.as_ref().unwrap();
                let bubbles_in_view = bbf.query(
                    self.view.chrom_idx,
                    self.view.view_start_bp as i32,
                    self.view.view_end_bp as i32,
                );
                let mut lane_ends: Vec<Vec<i32>> = vec![Vec::new(); n_samples];
                for b in bubbles_in_view {
                    let s = b.sample_idx as usize;
                    if s >= n_samples { continue; }
                    total_per_sample[s] += 1;
                    if !self.view.passes_filters(b) { continue; }
                    visible_per_sample[s] += 1;
                    total_visible_for_perf += 1;
                    let lane = lane_ends[s].iter().position(|&e| e < b.start);
                    if let Some(li) = lane {
                        lane_ends[s][li] = b.end;
                    } else {
                        lane_ends[s].push(b.end);
                    }
                }
                for (s, ends) in lane_ends.iter().enumerate() {
                    lanes_used_per_sample[s] = ends.len();
                }
            }
            self.last_query_ms = query_t0.elapsed().as_secs_f32() * 1000.0;
            self.last_visible_n = total_visible_for_perf;

            // Per-sample dynamic track height
            const SAMPLE_BASELINE_PAD: f32 = 22.0;
            let sample_heights: Vec<f32> = lanes_used_per_sample.iter()
                .map(|&n| {
                    let lanes = (n as f32).max(2.0); // always at least 2 visible
                    SAMPLE_BASELINE_PAD + (lanes + 1.0) * (LANE_H + LANE_GAP)
                })
                .collect();

            let avail_w = ui.available_width();

            // Per-sample visible height: clamp to max_lanes_per_sample.
            // The track's full content height is sample_heights[s]; the
            // visible portion (with mini-scrollbar) is visible_heights[s].
            let max_lanes_setting = self.view.max_lanes_per_sample as usize;
            let visible_heights: Vec<f32> = lanes_used_per_sample.iter()
                .enumerate()
                .map(|(s, &n)| {
                    let lanes_show = if max_lanes_setting == 0 { n.max(2) }
                                     else { n.min(max_lanes_setting).max(2) };
                    SAMPLE_BASELINE_PAD + (lanes_show as f32 + 1.0) * (LANE_H + LANE_GAP)
                        .min(sample_heights[s])
                })
                .collect();

            // No outer ScrollArea anymore — top tracks are fixed, each
            // sample-track has its own internal ScrollArea.
            //
            // Layout: ruler / seq / annot / sample[0] / sample[1] / ...
            // Pan/zoom interactions happen via raw pointer (ctx.input);
            // each individual response only contributes to hit-test.

            // Allocate the top header row (ruler + seq + annot).
            let n_annot = self.atlas.annot_tracks.len() as f32;
            let header_h = RULER_H + SEQ_H + n_annot * (ANNOT_H + 1.0);
            let (header_rect, header_response) = ui.allocate_exact_size(
                egui::Vec2::new(avail_w, header_h),
                egui::Sense::click_and_drag());
            // Cursor feedback: grab for the whole canvas-like region.
            if header_response.is_pointer_button_down_on() || header_response.dragged() {
                ctx.set_cursor_icon(egui::CursorIcon::Grabbing);
            } else if header_response.hovered() {
                ctx.set_cursor_icon(egui::CursorIcon::Grab);
            }
            // Use the header response as the canonical "canvas response"
            // for pan/zoom & for ruler-coordinate conversion.
            let response = header_response;
            let canvas_rect = header_rect;
            let p = ui.painter_at(canvas_rect);
            p.rect_filled(canvas_rect, 0.0, egui::Color32::from_gray(12));

                let bbf = self.atlas.bbf.as_ref().unwrap();
                let chrom_len_bp = bbf.chrom_index.get(self.view.chrom_idx as usize)
                    .map(|c| c.length_bp as f64);

                // ---- Vertical layout (inside ScrollArea) ----
                let mut y = canvas_rect.top();
                let track_left = canvas_rect.left() + SAMPLE_LABEL_W;
                let track_rect_w = canvas_rect.width() - SAMPLE_LABEL_W;

                // Ruler
                let ruler_rect = egui::Rect::from_min_size(
                    egui::Pos2::new(track_left, y),
                    egui::Vec2::new(track_rect_w, RULER_H));
                draw_ruler(&p, ruler_rect, &self.view);
                y += RULER_H;

                // Sequence track (always visible — placeholder when no FASTA)
                let seq_rect = egui::Rect::from_min_size(
                    egui::Pos2::new(track_left, y),
                    egui::Vec2::new(track_rect_w, SEQ_H));
                let chrom_name = bbf.chroms.strings.get(self.view.chrom_idx as usize)
                    .cloned().unwrap_or_default();
                if let Some(fa) = &self.atlas.reference {
                    draw_sequence(&p, seq_rect, &self.view, fa, &chrom_name);
                } else {
                    p.rect_filled(seq_rect, 0.0, egui::Color32::from_gray(20));
                    p.text(
                        seq_rect.center(),
                        egui::Align2::CENTER_CENTER,
                        "no reference loaded — Open Reference (.fa) to see ATCG at deep zoom",
                        egui::FontId::proportional(10.0),
                        egui::Color32::from_gray(140),
                    );
                }
                p.text(
                    egui::Pos2::new(canvas_rect.left() + 4.0, seq_rect.center().y),
                    egui::Align2::LEFT_CENTER,
                    "sequence",
                    egui::FontId::monospace(9.0),
                    egui::Color32::from_gray(160));
                y += SEQ_H;

                // Annotation tracks
                for at in &self.atlas.annot_tracks {
                    let annot_rect = egui::Rect::from_min_size(
                        egui::Pos2::new(track_left, y),
                        egui::Vec2::new(track_rect_w, ANNOT_H));
                    draw_annotation(&p, annot_rect, &self.view, &chrom_name, at);
                    p.text(
                        egui::Pos2::new(canvas_rect.left() + 4.0, annot_rect.center().y),
                        egui::Align2::LEFT_CENTER,
                        &at.label,
                        egui::FontId::proportional(11.0),
                        egui::Color32::from_gray(160));
                    y += ANNOT_H + 1.0;
                }

                // Bubble tracks: one per sample, each in its OWN ScrollArea.
                let bp_per_px = self.view.bp_per_px(track_rect_w);
                let auto_mode = pick_mode(bp_per_px);
                let mode = if self.view.force_haplotype_mode {
                    crate::render::BubbleMode::Haplotype
                } else { auto_mode };
                let mut new_hover = None::<usize>;
                let mut pending_msa_idx: Option<usize> = None;
                let render_t0 = std::time::Instant::now();

            for s in 0..n_samples {
                let inner_h = sample_heights[s];
                let visible_h = visible_heights[s].min(inner_h);
                ui.horizontal(|ui| {
                    // -- left gutter (fixed width, shows sample label / counts) --
                    let (gutter_rect, _) = ui.allocate_exact_size(
                        egui::Vec2::new(SAMPLE_LABEL_W, visible_h),
                        egui::Sense::hover());
                    let gp = ui.painter_at(gutter_rect);
                    gp.rect_filled(gutter_rect, 0.0, egui::Color32::from_gray(14));
                    let bbf_ref = self.atlas.bbf.as_ref().unwrap();
                    let label = bbf_ref.samples.get(s as u32).unwrap_or("?").to_string();
                    gp.text(
                        egui::Pos2::new(gutter_rect.left() + 6.0, gutter_rect.bottom() - 36.0),
                        egui::Align2::LEFT_CENTER,
                        &label,
                        egui::FontId::proportional(12.0),
                        egui::Color32::from_gray(220));
                    gp.text(
                        egui::Pos2::new(gutter_rect.left() + 6.0, gutter_rect.bottom() - 22.0),
                        egui::Align2::LEFT_CENTER,
                        format!("{:?} | lanes={}", mode, lanes_used_per_sample[s]),
                        egui::FontId::monospace(9.0),
                        egui::Color32::from_gray(120));
                    gp.text(
                        egui::Pos2::new(gutter_rect.left() + 6.0, gutter_rect.bottom() - 10.0),
                        egui::Align2::LEFT_CENTER,
                        format!("{}/{}", visible_per_sample[s], total_per_sample[s]),
                        egui::FontId::monospace(10.0),
                        egui::Color32::from_rgb(0xfd, 0xae, 0x61));

                    // -- right: ScrollArea with the bubble track --
                    egui::ScrollArea::vertical()
                        .max_height(visible_h)
                        .id_source(format!("sample_track_{}", s))
                        .auto_shrink([false; 2])
                        .show(ui, |ui| {
                            let avail_w_inner = ui.available_width();
                            let (track_rect, track_resp) = ui.allocate_exact_size(
                                egui::Vec2::new(avail_w_inner, inner_h),
                                egui::Sense::click_and_drag());
                            let pp = ui.painter_at(track_rect);
                            pp.rect_filled(track_rect, 0.0, egui::Color32::from_gray(12));
                            let local_hover = track_resp.hover_pos();
                            let bbf_inner = self.atlas.bbf.as_ref().unwrap();
                            let hit = draw_bubble_track(&pp, track_rect, &self.view,
                                &self.atlas, bbf_inner, s as u8, local_hover, mode,
                                self.view.selected_bubble);
                            if let Some(h) = hit {
                                if new_hover.is_none() {
                                    new_hover = Some(h.global_index);
                                }
                            }
                            // also accept right-click and clicks here — but we
                            // route them up to the central state via the same
                            // hover_bubble / selected_bubble mechanism.
                            if track_resp.clicked() && new_hover.is_some() {
                                self.view.selected_bubble = new_hover;
                            }
                            // right-click context menu attached to this track
                            track_resp.context_menu(|ui| {
                                if let Some(idx) = new_hover {
                                    if let Some(rec) = bbf_inner.bubbles.get(idx) {
                                        let bname = bbf_inner.bubble_names.get(rec.bubble_name_idx)
                                            .unwrap_or("(no name)").to_string();
                                        let sample = bbf_inner.samples.get(rec.sample_idx as u32)
                                            .unwrap_or("?").to_string();
                                        let class = bbf_inner.classes.get(rec.class_idx as u32)
                                            .unwrap_or("?").to_string();
                                        let chrom = bbf_inner.chroms.get(rec.chrom_idx as u32)
                                            .unwrap_or("?").to_string();
                                        let dbvar_id = bbf_inner.dbvar_ids.get(rec.dbvar_id_idx)
                                            .unwrap_or("").to_string();
                                        ui.label(egui::RichText::new(format!("{}  ({})", &bname, &sample)).strong());
                                        ui.label(format!("{}:{}..{}  len={}  VAF={:.4}",
                                            &chrom, rec.start, rec.end, rec.length(), rec.vaf()));
                                        ui.label(format!("class: {}{}", &class,
                                            if rec.is_shared() { "  (shared ≥2 samples)" } else { "" }));
                                        let n_alts = rec.n_alts as u32;
                                        let total_reads = rec.total_reads_approx();
                                        ui.label(format!(
                                            "read support: ≈{} reads across {} alts → minor alt ≈{} reads ({:.2}% VAF)",
                                            total_reads, n_alts,
                                            (rec.vaf() * total_reads as f32).round() as u32,
                                            rec.vaf() * 100.0));
                                        if !dbvar_id.is_empty() {
                                            ui.label(format!("dbVar match: {}  reciprocal-overlap={:.2}",
                                                &dbvar_id, rec.dbvar_recip()));
                                        }
                                        if rec.pop_count() > 0 {
                                            let interp = match rec.pop_count() {
                                                1        => "singleton in dbVar",
                                                2..=5    => "rare (≤5 nssv obs.)",
                                                6..=50   => "uncommon",
                                                51..=500 => "common",
                                                _        => "very common",
                                            };
                                            ui.label(format!(
                                                "dbVar pop count: {} nssv  ({})",
                                                rec.pop_count(), interp));
                                        } else {
                                            ui.label("dbVar pop count: 0  (novel — no region match)");
                                        }
                                        // .bbz alt-sequences (if loaded for this bubble)
                                        let alts_for_bubble = self.atlas.bbz.as_ref()
                                            .and_then(|z| z.alts_for(&bname));
                                        if let Some(alts) = alts_for_bubble {
                                            ui.label(format!("alt-sequences in .bbz: {} alts", alts.len()));
                                        }
                                        ui.separator();
                                        if ui.button("Copy bubble name").clicked() {
                                            ui.output_mut(|o| o.copied_text = bname.clone());
                                            ui.close_menu();
                                        }
                                        if ui.button("Copy as BED line").clicked() {
                                            let bed = format!("{}\t{}\t{}\t{}\t{:.4}",
                                                &chrom, rec.start, rec.end, &bname, rec.vaf());
                                            ui.output_mut(|o| o.copied_text = bed);
                                            ui.close_menu();
                                        }
                                        if ui.button("Copy as TSV row").clicked() {
                                            let tsv = format!("{sample}\t{bname}\t{chrom}\t{s2}\t{e}\t{vaf:.4}\t{class}\t{dbid}\t{recip:.3}",
                                                s2 = rec.start, e = rec.end,
                                                vaf = rec.vaf(), recip = rec.dbvar_recip(),
                                                dbid = &dbvar_id);
                                            ui.output_mut(|o| o.copied_text = tsv);
                                            ui.close_menu();
                                        }
                                        ui.separator();
                                        // External lookups: dbVar (if matched) + UCSC by coords.
                                        // nssv accessions go straight to the variant page; bare
                                        // region IDs (nstd) need the search endpoint. If no
                                        // dbVar id, fall back to a coordinate search.
                                        let dbvar_url = if !dbvar_id.is_empty() {
                                            if dbvar_id.starts_with("nssv") {
                                                format!("https://www.ncbi.nlm.nih.gov/dbvar/variants/{}/", dbvar_id)
                                            } else {
                                                format!("https://www.ncbi.nlm.nih.gov/dbvar/?term={}", dbvar_id)
                                            }
                                        } else {
                                            format!("https://www.ncbi.nlm.nih.gov/dbvar/?term=({}%3A{}-{})%5BBase+Position%5D",
                                                &chrom, rec.start, rec.end)
                                        };
                                        let dbvar_label = if dbvar_id.is_empty() {
                                            "Search dbVar by coordinates ↗".to_string()
                                        } else {
                                            format!("Open {} in dbVar ↗", &dbvar_id)
                                        };
                                        if ui.button(dbvar_label).clicked() {
                                            ui.ctx().output_mut(|o| {
                                                o.open_url = Some(egui::output::OpenUrl {
                                                    url: dbvar_url,
                                                    new_tab: true,
                                                });
                                            });
                                            ui.close_menu();
                                        }
                                        // UCSC quick-jump by coordinates (defaults to GRCh38).
                                        if ui.button("View locus in UCSC (hg38) ↗").clicked() {
                                            let ucsc = format!(
                                                "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={}%3A{}-{}",
                                                &chrom, rec.start.max(0), rec.end.max(rec.start + 1));
                                            ui.ctx().output_mut(|o| {
                                                o.open_url = Some(egui::output::OpenUrl {
                                                    url: ucsc,
                                                    new_tab: true,
                                                });
                                            });
                                            ui.close_menu();
                                        }
                                        ui.separator();
                                        // .bbz-aware actions
                                        let bbz_loaded = self.atlas.bbz.is_some();
                                        let alts_avail = alts_for_bubble.map(|a| a.len()).unwrap_or(0);
                                        if ui.add_enabled(bbz_loaded && alts_avail > 0,
                                            egui::Button::new(format!("Export {} ALTs as FASTA…", alts_avail))).clicked()
                                        {
                                            if let Some(alts) = alts_for_bubble {
                                                let mut fasta = String::new();
                                                for (i, a) in alts.iter().enumerate() {
                                                    let path_str: Vec<String> = a.path_nodes.iter()
                                                        .map(|n| n.to_string()).collect();
                                                    fasta.push_str(&format!(
                                                        ">{bname}_alt{i} sample={sample} {chrom}:{s2}-{e} len={l} path={p} \n",
                                                        i = i, s2 = rec.start, e = rec.end,
                                                        l = a.seq.len(), p = path_str.join(",")));
                                                    // wrap to 60 cols
                                                    for chunk in a.seq.chunks(60) {
                                                        fasta.push_str(std::str::from_utf8(chunk).unwrap_or(""));
                                                        fasta.push('\n');
                                                    }
                                                }
                                                if let Some(p) = rfd::FileDialog::new()
                                                    .add_filter("FASTA", &["fa","fasta"])
                                                    .set_file_name(format!("{}.fa", &bname))
                                                    .save_file()
                                                {
                                                    let _ = std::fs::write(&p, fasta);
                                                }
                                            }
                                            ui.close_menu();
                                        }
                                        if ui.add_enabled(bbz_loaded && alts_avail > 0,
                                            egui::Button::new("Copy first ALT sequence to clipboard")).clicked()
                                        {
                                            if let Some(alts) = alts_for_bubble {
                                                if let Some(a0) = alts.first() {
                                                    let s = std::str::from_utf8(&a0.seq).unwrap_or("").to_string();
                                                    ui.output_mut(|o| o.copied_text = s);
                                                }
                                            }
                                            ui.close_menu();
                                        }
                                        if ui.add_enabled(bbz_loaded,
                                            egui::Button::new("Show MSA at this locus…")).clicked()
                                        {
                                            pending_msa_idx = Some(idx);
                                            ui.close_menu();
                                        }
                                    } else {
                                        ui.label("hover over a bubble first");
                                    }
                                }
                            });
                        });
                });
            }
            self.view.hover_bubble = new_hover;
            self.last_render_ms = render_t0.elapsed().as_secs_f32() * 1000.0;
            // ---- pending MSA open (set from a context-menu click) ----
            if let Some(idx) = pending_msa_idx {
                self.open_msa_at_locus(idx);
            }

            // ---- variables consumed by the trailing interaction / overlay code ----
            let bbf = self.atlas.bbf.as_ref().unwrap();
            let track_left = canvas_rect.left() + SAMPLE_LABEL_W;
            let ruler_rect = egui::Rect::from_min_size(
                egui::Pos2::new(track_left, canvas_rect.top()),
                egui::Vec2::new(canvas_rect.width() - SAMPLE_LABEL_W, RULER_H));
            let track_rect_w = ruler_rect.width();

            // Unmapped overlay drawn on the *header* rect only (ruler / seq /
            // annot).  Per-sample tracks live in their own ScrollAreas and
            // the overlay covering them gets messy; the header overlay alone
            // is enough to communicate "this region is past the analysis".
            if let Some((a_start, a_end)) = self.atlas.analysis_range(self.view.chrom_idx) {
                let xs = bp_to_x(a_start as f64, &self.view, ruler_rect);
                let xe = bp_to_x(a_end as f64, &self.view, ruler_rect);
                if xs > track_left {
                    let r = egui::Rect::from_min_max(
                        egui::Pos2::new(track_left, canvas_rect.top()),
                        egui::Pos2::new(xs.min(canvas_rect.right()), canvas_rect.bottom()));
                    if r.width() > 0.5 {
                        draw_unmapped_overlay(&p, r,
                            &format!("unmapped (before {})", format_bp(a_start as i64)));
                    }
                }
                if xe < canvas_rect.right() {
                    let r = egui::Rect::from_min_max(
                        egui::Pos2::new(xe.max(track_left), canvas_rect.top()),
                        egui::Pos2::new(canvas_rect.right(), canvas_rect.bottom()));
                    if r.width() > 0.5 {
                        draw_unmapped_overlay(&p, r,
                            &format!("unmapped (after {})", format_bp(a_end as i64)));
                    }
                }
            }

            // ----- interaction: pan via raw pointer tracking -----
            // We bypass response.dragged() because some Linux/Wayland setups
            // don't fire it reliably with allocate_exact_size, AND so we can
            // explicitly suppress scroll-zoom while the user is panning.
            let primary_down = ctx.input(|i| i.pointer.primary_down());
            let pointer_delta = ctx.input(|i| i.pointer.delta());
            let pointer_pos = ctx.input(|i| i.pointer.hover_pos());

            // Track our own drag state: drag starts when primary goes down on
            // the canvas, ends when button is released.
            if primary_down && response.is_pointer_button_down_on() {
                self.drag_active = true;
            }
            if !primary_down {
                self.drag_active = false;
            }

            if self.drag_active && pointer_delta.x.abs() > 0.0 {
                let bp_delta = -(pointer_delta.x as f64) * bp_per_px;
                self.view.pan(bp_delta);
                self.view.clamp(chrom_len_bp);
            }

            // Zoom requires Ctrl/Cmd modifier so that plain wheel scrolling
            // inside a per-sample ScrollArea continues to scroll lanes
            // without simultaneously zooming the genomic view.
            let zoom_modifier = ctx.input(|i| i.modifiers.ctrl || i.modifiers.command);
            if !self.drag_active && zoom_modifier && scroll_y.abs() > 0.0 {
                let hp = pointer_pos.or_else(|| response.hover_pos());
                if let Some(hp) = hp {
                    let pivot_bp = x_to_bp(hp.x, &self.view, ruler_rect);
                    let factor = if scroll_y > 0.0 { 0.85 } else { 1.0 / 0.85 };
                    self.view.zoom(factor, pivot_bp);
                    self.view.clamp(chrom_len_bp);
                }
            }
            if key_plus {
                let pivot = (self.view.view_start_bp + self.view.view_end_bp) / 2.0;
                self.view.zoom(0.7, pivot);
                self.view.clamp(chrom_len_bp);
            }
            if key_minus {
                let pivot = (self.view.view_start_bp + self.view.view_end_bp) / 2.0;
                self.view.zoom(1.0 / 0.7, pivot);
                self.view.clamp(chrom_len_bp);
            }
            if key_left {
                self.view.pan(-self.view.span() * 0.2);
                self.view.clamp(chrom_len_bp);
            }
            if key_right {
                self.view.pan(self.view.span() * 0.2);
                self.view.clamp(chrom_len_bp);
            }
            if key_f {
                if let Some(len) = chrom_len_bp {
                    self.view.view_start_bp = 0.0;
                    self.view.view_end_bp = len;
                }
            }

            // selected bubble highlight info — just a small overlay
            if let Some(idx) = self.view.selected_bubble {
                if let Some(rec) = bbf.bubbles.get(idx) {
                    let info = format!(
                        "Selected: {} {} VAF={:.4} length={} dbVar={} recip={:.2}",
                        bbf.samples.get(rec.sample_idx as u32).unwrap_or("?"),
                        bbf.bubble_names.get(rec.bubble_name_idx).unwrap_or(""),
                        rec.vaf(),
                        format_bp(rec.length() as i64),
                        bbf.dbvar_ids.get(rec.dbvar_id_idx).unwrap_or(""),
                        rec.dbvar_recip(),
                    );
                    p.text(
                        egui::Pos2::new(canvas_rect.left() + 4.0, canvas_rect.bottom() - 14.0),
                        egui::Align2::LEFT_BOTTOM,
                        info,
                        egui::FontId::monospace(11.0),
                        egui::Color32::YELLOW,
                    );
                }
            }

            // ----- Perf overlay (toggle with F11 or --perf flag) -----
            if self.perf_visible {
                let txt = format!(
                    "frame={:.1}ms  query={:.2}ms  render={:.1}ms  visible={}  bp/px={:.1}",
                    self.last_frame_ms, self.last_query_ms,
                    self.last_render_ms, self.last_visible_n,
                    bp_per_px,
                );
                let pos = egui::Pos2::new(canvas_rect.right() - 6.0, canvas_rect.top() + 4.0);
                p.rect_filled(
                    egui::Rect::from_min_size(
                        egui::Pos2::new(pos.x - 380.0, pos.y - 1.0),
                        egui::Vec2::new(380.0, 16.0)),
                    2.0,
                    egui::Color32::from_rgba_premultiplied(0, 0, 0, 200));
                p.text(pos, egui::Align2::RIGHT_TOP, txt,
                    egui::FontId::monospace(10.0),
                    egui::Color32::from_rgb(0x6f, 0xe0, 0xff));
            }
        });

        // ---- MSA pop-up window (if open) ----
        if self.msa_open {
            if let Some(state) = self.msa_state.as_mut() {
                let mut open = self.msa_open;
                msa::render_msa_window(ctx, state, &mut open);
                self.msa_open = open;
            }
        }

        self.last_frame_ms = frame_t0.elapsed().as_secs_f32() * 1000.0;
        if self.perf_visible {
            ctx.request_repaint();   // keep redrawing so frame count updates
        }
    }
}

fn main() -> eframe::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let opts = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1400.0, 800.0])
            .with_min_inner_size([900.0, 500.0])
            .with_title("VariantPaths"),
        ..Default::default()
    };
    eframe::run_native(
        "VariantPaths",
        opts,
        Box::new(move |cc| Box::new(App::new(cc, args.clone()))),
    )
}
