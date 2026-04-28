//! Viewport state: current chromosome, zoom, filters.

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum SharedFilter {
    All,
    SharedOnly,
    PrivateOnly,
}

#[derive(Clone)]
pub struct ViewState {
    pub chrom_idx: u16,
    pub view_start_bp: f64,
    pub view_end_bp: f64,

    // ----- user-controllable filters -----
    pub sample_visible: Vec<bool>,
    pub class_visible: Vec<bool>,
    pub vaf_min: f32,            // inclusive, default 0.0
    pub vaf_max: f32,            // inclusive, default 1.0
    pub min_total_reads_log: u8, // log2(R)*8 quantized, default 0 (no filter)
    pub min_length_bp: i32,      // default 0
    pub max_length_bp: i32,      // default i32::MAX
    pub shared_filter: SharedFilter,
    /// Max lanes shown per sample-track before its inner scrollbar kicks
    /// in.  0 = unlimited (show all lanes, no per-sample scroll).
    pub max_lanes_per_sample: u32,

    pub jump_input: String,
    pub selected_bubble: Option<usize>,
    pub hover_bubble: Option<usize>,
}

impl ViewState {
    pub fn fit_chrom(chrom_length: u32) -> Self {
        Self {
            chrom_idx: 0,
            view_start_bp: 0.0,
            view_end_bp: chrom_length.max(1) as f64,
            sample_visible: vec![],
            class_visible: vec![],
            vaf_min: 0.0,
            vaf_max: 1.0,
            min_total_reads_log: 0,
            min_length_bp: 0,
            max_length_bp: i32::MAX,
            shared_filter: SharedFilter::All,
            max_lanes_per_sample: 6,
            jump_input: String::new(),
            selected_bubble: None,
            hover_bubble: None,
        }
    }

    /// Returns true if this bubble passes all user-set filters.
    pub fn passes_filters(&self, b: &crate::format::bbf::BubbleRec) -> bool {
        // sample visibility (vec may be shorter than n_samples right after load)
        if let Some(&v) = self.sample_visible.get(b.sample_idx as usize) {
            if !v { return false; }
        }
        if let Some(&v) = self.class_visible.get(b.class_idx as usize) {
            if !v { return false; }
        }
        let v = b.vaf();
        if v < self.vaf_min || v > self.vaf_max { return false; }
        if b.total_reads_log < self.min_total_reads_log { return false; }
        let len = b.length();
        if len < self.min_length_bp { return false; }
        if len > self.max_length_bp { return false; }
        match self.shared_filter {
            SharedFilter::All => {}
            SharedFilter::SharedOnly => if !b.is_shared() { return false; }
            SharedFilter::PrivateOnly => if b.is_shared() { return false; }
        }
        true
    }

    pub fn span(&self) -> f64 {
        (self.view_end_bp - self.view_start_bp).max(1.0)
    }

    /// Zoom by `factor` around `pivot_bp`.  factor < 1 zoom in, > 1 zoom out.
    pub fn zoom(&mut self, factor: f64, pivot_bp: f64) {
        let new_span = (self.span() * factor).max(10.0); // can't zoom past 10 bp
        let pivot_frac = (pivot_bp - self.view_start_bp) / self.span();
        self.view_start_bp = pivot_bp - pivot_frac * new_span;
        self.view_end_bp = self.view_start_bp + new_span;
        self.clamp(None);
    }

    /// Pan in bp.
    pub fn pan(&mut self, delta_bp: f64) {
        self.view_start_bp += delta_bp;
        self.view_end_bp += delta_bp;
        self.clamp(None);
    }

    /// Soft clamp.  We let the user pan past the analyzed region — the
    /// renderer will paint an "unmapped" overlay there.  We only protect
    /// against degenerate spans and runaway offsets.
    pub fn clamp(&mut self, _chrom_len_bp: Option<f64>) {
        // Never let span collapse to <= 0.
        if self.view_end_bp <= self.view_start_bp {
            self.view_end_bp = self.view_start_bp + 1.0;
        }
        // Generous safety cap so floating-point doesn't go nuts.
        const MAX_OFFSET: f64 = 5.0e9;
        if self.view_start_bp < -MAX_OFFSET {
            let s = self.span();
            self.view_start_bp = -MAX_OFFSET;
            self.view_end_bp = self.view_start_bp + s;
        }
        if self.view_end_bp > MAX_OFFSET {
            let s = self.span();
            self.view_end_bp = MAX_OFFSET;
            self.view_start_bp = self.view_end_bp - s;
        }
    }

    pub fn bp_per_px(&self, canvas_width_px: f32) -> f64 {
        self.span() / (canvas_width_px.max(1.0) as f64)
    }
}
