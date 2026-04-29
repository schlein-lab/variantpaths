//! Top toolbar: file open buttons, chrom switcher, jump-to, filter sliders.
//!
//! Extracted from `main.rs::App::update` so the file/nav row and the rather
//! sprawling filter panel don't drown out the actual rendering pipeline in
//! the main update loop.  Behaviour is unchanged.

use crate::App;
use crate::render::format_bp;
use crate::view::SharedFilter;

pub fn show_topbar(app: &mut App, ctx: &egui::Context) {
    egui::TopBottomPanel::top("topbar").show(ctx, |ui| {
        // ----- First row: file open + nav -----
        ui.horizontal(|ui| {
            if ui.button("Open .vpf").clicked() {
                if let Some(p) = rfd::FileDialog::new()
                    .add_filter("VariantPaths File", &["vpf"])
                    .pick_file() { app.load_bbf(&p); }
            }
            if ui.button("Open Reference (.fa)").clicked() {
                if let Some(p) = rfd::FileDialog::new()
                    .add_filter("FASTA", &["fa", "fasta", "fna"])
                    .pick_file() { app.load_fasta(&p); }
            }
            if ui.button("Open BED").clicked() {
                if let Some(p) = rfd::FileDialog::new()
                    .add_filter("BED", &["bed"])
                    .pick_file() { app.load_bed(&p); }
            }
            if ui.button("Open .vpz (sequences)").clicked() {
                if let Some(p) = rfd::FileDialog::new()
                    .add_filter("VariantPaths Sequences", &["vpz"])
                    .pick_file() { app.load_bbz(&p); }
            }
            ui.separator();
            if ui.add_enabled(
                app.atlas.bbf.is_some(),
                egui::Button::new("Export filtered TSV…"),
            ).clicked() {
                app.export_filtered_tsv();
            }
            ui.separator();
            ui.checkbox(&mut app.heatmap_visible, "Heatmap");
            ui.separator();

            let has_bbf = app.atlas.bbf.is_some();
            if has_bbf {
                let (chrom_names, chrom_lens, current_idx) = {
                    let b = app.atlas.bbf.as_ref().unwrap();
                    (
                        b.chroms.strings.clone(),
                        b.chrom_index.iter().map(|c| c.length_bp).collect::<Vec<_>>(),
                        app.view.chrom_idx,
                    )
                };
                let mut new_chrom_idx: Option<u16> = None;
                egui::ComboBox::from_label("chrom")
                    .selected_text(chrom_names.get(current_idx as usize).cloned().unwrap_or_default())
                    .show_ui(ui, |ui| {
                        for (i, name) in chrom_names.iter().enumerate() {
                            if ui.selectable_label(i as u16 == current_idx, name).clicked() {
                                new_chrom_idx = Some(i as u16);
                            }
                        }
                    });
                if let Some(i) = new_chrom_idx {
                    app.view.chrom_idx = i;
                    let len = chrom_lens.get(i as usize).copied().unwrap_or(1);
                    app.view.view_start_bp = 0.0;
                    app.view.view_end_bp = len as f64;
                }
                ui.add_sized(
                    [240.0, 20.0],
                    egui::TextEdit::singleline(&mut app.view.jump_input)
                        .hint_text("chr14:105_700_000  ·  IGHG4  ·  chr14"),
                );
                if ui.button("Go").clicked() {
                    let s = app.view.jump_input.clone();
                    app.jump_to(&s);
                }
                ui.label(format!(
                    "{}:{} - {}",
                    chrom_names.get(app.view.chrom_idx as usize).map(|s| s.as_str()).unwrap_or(""),
                    format_bp(app.view.view_start_bp as i64),
                    format_bp(app.view.view_end_bp as i64),
                ));
            } else {
                ui.label("no .vpf loaded — drop a file or use Open .vpf");
            }
        });

        // ----- Second row: filters (collapsible) -----
        let Some(b) = app.atlas.bbf.as_ref() else { return; };
        let sample_names = b.samples.strings.clone();
        let class_names = b.classes.strings.clone();
        let _ = b;

        ui.horizontal_wrapped(|ui| {
            egui::CollapsingHeader::new("Filters")
                .default_open(true)
                .show(ui, |ui| {
                ui.horizontal_wrapped(|ui| {
                    ui.label("VAF:");
                    ui.add(egui::Slider::new(&mut app.view.vaf_min, 0.0..=1.0)
                        .text("min").logarithmic(true).min_decimals(4).step_by(0.0001));
                    ui.add(egui::Slider::new(&mut app.view.vaf_max, 0.0..=1.0)
                        .text("max").logarithmic(true).min_decimals(4).step_by(0.0001));
                    if app.view.vaf_min > app.view.vaf_max {
                        app.view.vaf_max = app.view.vaf_min;
                    }
                    if ui.small_button("subclonal <5%").clicked() {
                        app.view.vaf_min = 0.0; app.view.vaf_max = 0.05;
                    }
                    if ui.small_button("intermediate 5-25%").clicked() {
                        app.view.vaf_min = 0.05; app.view.vaf_max = 0.25;
                    }
                    if ui.small_button("germline-like ≥25%").clicked() {
                        app.view.vaf_min = 0.25; app.view.vaf_max = 1.0;
                    }
                    if ui.small_button("all").clicked() {
                        app.view.vaf_min = 0.0; app.view.vaf_max = 1.0;
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("min total_reads:");
                    let mut approx_reads = if app.view.min_total_reads_log == 0 { 0u32 }
                        else { (2.0_f32).powf(app.view.min_total_reads_log as f32 / 8.0) as u32 };
                    if ui.add(egui::Slider::new(&mut approx_reads, 0..=500).text("reads"))
                        .changed()
                    {
                        app.view.min_total_reads_log = if approx_reads == 0 { 0 }
                            else { ((approx_reads as f32).log2() * 8.0).round().clamp(0.0, 255.0) as u8 };
                    }
                    ui.separator();
                    ui.label("length range (bp):");
                    ui.add(egui::DragValue::new(&mut app.view.min_length_bp)
                        .speed(100).clamp_range(0..=1_000_000_000));
                    ui.label("–");
                    ui.add(egui::DragValue::new(&mut app.view.max_length_bp)
                        .speed(100).clamp_range(0..=i32::MAX));
                    if ui.small_button("reset").clicked() {
                        app.view.min_length_bp = 0;
                        app.view.max_length_bp = i32::MAX;
                        app.view.min_total_reads_log = 0;
                    }
                });

                if app.view.sample_visible.len() != sample_names.len() {
                    app.view.sample_visible = vec![true; sample_names.len()];
                }
                if app.view.class_visible.len() != class_names.len() {
                    app.view.class_visible = vec![true; class_names.len()];
                }

                ui.horizontal_wrapped(|ui| {
                    ui.label("recurrence:");
                    ui.radio_value(&mut app.view.shared_filter, SharedFilter::All, "all");
                    ui.radio_value(&mut app.view.shared_filter, SharedFilter::SharedOnly,
                        "shared (≥2 samples)");
                    ui.radio_value(&mut app.view.shared_filter, SharedFilter::PrivateOnly,
                        "private (1 sample)");
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("max lanes per sample:");
                    let mut v = app.view.max_lanes_per_sample;
                    ui.add(egui::Slider::new(&mut v, 0u32..=50).integer().text("(0 = all)"));
                    app.view.max_lanes_per_sample = v;
                    if ui.small_button("6").clicked() { app.view.max_lanes_per_sample = 6; }
                    if ui.small_button("12").clicked() { app.view.max_lanes_per_sample = 12; }
                    if ui.small_button("all").clicked() { app.view.max_lanes_per_sample = 0; }
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("max dbVar pop count:");
                    let mut cap_int = if app.view.max_pop_count == u16::MAX { 10000u32 }
                        else { app.view.max_pop_count as u32 };
                    let resp = ui.add(egui::Slider::new(&mut cap_int, 0u32..=10000)
                        .integer().logarithmic(true).text("nssv"));
                    if resp.changed() {
                        app.view.max_pop_count = if cap_int >= 10000 { u16::MAX }
                            else { cap_int as u16 };
                    }
                    if ui.small_button("rare ≤5").clicked() { app.view.max_pop_count = 5; }
                    if ui.small_button("≤50").clicked() { app.view.max_pop_count = 50; }
                    if ui.small_button("∞").clicked() { app.view.max_pop_count = u16::MAX; }
                });
                ui.horizontal_wrapped(|ui| {
                    ui.checkbox(&mut app.view.force_haplotype_mode,
                        "haplotype mode (two parallel lines, bubbles routed by VAF≷0.5)");
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("samples:");
                    for (i, n) in sample_names.iter().enumerate() {
                        ui.checkbox(&mut app.view.sample_visible[i], n);
                    }
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("classes:");
                    for (i, n) in class_names.iter().enumerate() {
                        let col = app.atlas.class_color(i as u8);
                        let label = egui::RichText::new(n).color(col);
                        ui.checkbox(&mut app.view.class_visible[i], label);
                    }
                });
            });
        });
    });
}
