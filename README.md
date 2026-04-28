# VariantPaths

**Standalone genome-graph viewer for structural variant paths.**
A single ~8 MB executable that reads a portable `.vpf` data file and renders
genome-wide structural variants as graph paths — stacked like reads in IGV,
with optional alt-sequences, multiple-sequence alignment, and reference
overlay. No installer. No server. No browser plugin.

🌐 [variantpaths.com](https://variantpaths.com)

---

## Why VariantPaths

Existing genome browsers (IGV, JBrowse, UCSC) show reads against a single
linear reference. Once you have **graph-based assembly output** — bubbles,
alt-paths, paralog rearrangements — those tools fall short:

- A bubble with two alt-paths doesn't reduce to a single VCF line
- Paralog-rich loci (IGH, KIR, MHC) stack tens of overlapping variants
- Per-cell mosaic CNV at <5 % VAF is invisible in CIGAR-string views

**VariantPaths** treats every alt-path as a first-class object: lay them out
as horizontal **lanes** (IGV-style pile-up), annotate each with VAF + dbVar
match, optionally show the underlying ALT nucleotide sequence, and run an
on-the-spot multiple-sequence alignment of all alt-paths at a locus.

## Features

- **Browse**: chr-jump, mouse-pan, ctrl+scroll-zoom, vertical lane scroll per
  sample
- **Filter**: VAF range, min reads, length, classification, recurrence
  (private vs shared), max-lanes-per-sample
- **Inspect**: hover for metadata; right-click for clipboard / FASTA /
  TSV-row export
- **Sequence**: load a `.vpz` companion to see ALT sequences inline; load a
  reference FASTA to see ATCG at <5 bp/px
- **Compare**: "Show MSA at this locus" — pop-up window with all bubble
  alt-paths aligned (Needleman-Wunsch on demand), divergent columns
  highlighted, copy-to-FASTA
- **Annotate**: built-in IGH gene track, custom BED loadable
- **Cross-platform**: native Linux / Windows `.exe` (cross-compiled via
  cross-rs); ~8 MB statically linked, no runtime dependencies
- **Email-able data**: `.vpf` for 20 000 IGH bubbles ≈ 250 KB, `.vpz` with
  alt-sequences ≈ 75 KB

## File formats

| Extension | Content | Typical size |
|---|---|---|
| `.vpf` | bubble topology + classification + dbVar match (binary, zstd) | 250 KB / 20k bubbles |
| `.vpz` | alt-path nucleotide sequences (binary, zstd) | 75 KB / 1k alt-paths |

The `.vpf` is mandatory for the viewer; `.vpz` unlocks the sequence panel,
FASTA export, and MSA viewer; a reference FASTA unlocks the `ATCG` track at
deep zoom.

## Build

```sh
# requires: rustup (1.76+), python3 + zstandard (only for build_vpf.py)
git clone https://github.com/schlein-lab/variantpaths
cd variantpaths
cargo build --release
./target/release/variantpaths samples/igh_3sample.vpf samples/igh_3sample.vpz
```

### Cross-compile to Windows

```sh
sudo apt install -y mingw-w64    # one-time
rustup target add x86_64-pc-windows-gnu
./build_windows.sh               # uses cross-rs (Docker)
# → dist/VariantPaths.exe (~4 MB statically linked PE32+)
```

## Generate your own `.vpf` / `.vpz`

```sh
# from a TSV with bubble data (BRANCH atlas output works directly):
python3 build_vpf.py \
    --tsv your_atlas.tsv \
    --out your.vpf \
    --reference GRCh38

# from per-sample GFAs + the same TSV (alt-sequence extraction):
python3 build_vpz.py \
    --tsv your_atlas.tsv \
    --gfa-map "sampleA=path/to/A.gfa,sampleB=path/to/B.gfa" \
    --max-paths-per-bubble 20 \
    --max-alt-len 50000 \
    --out your.vpz
```

The `.vpf` schema lives in [`src/format/bbf.rs`](src/format/bbf.rs);
the `.vpz` schema in [`src/format/bbz.rs`](src/format/bbz.rs).
*(File-format type names retain `bbf`/`bbz` from the project's prototype
phase; the on-disk extensions are `.vpf`/`.vpz`.)*

## Keyboard / mouse

| Action | Input |
|---|---|
| Pan genomic position | left-mouse drag |
| Zoom in / out | **Ctrl + mouse wheel** (centered on cursor) |
| Vertical lane scroll | mouse wheel inside a sample track |
| Fit chromosome | `f` |
| Step pan | `←` / `→` |
| Zoom in / out (centered) | `+` / `-` |
| Toggle perf overlay | `F11` |
| Right-click bubble | context menu (copy, export FASTA, MSA…) |

## Status

This is an **active research codebase** developed by the Schlein Lab. The
data formats are versioned (`VPF1`, `VPZ1`); we will avoid breaking format
changes without a version bump. APIs in `src/` are not yet stable.

Originally built for IGHG4 single-cell CNV instability analysis on HPRC
LCL samples; designed to generalize to whole-genome graph SV analysis.

## License

[MIT](LICENSE)

## Citation

If you use VariantPaths in published work, please cite *(citation TBD —
preprint forthcoming)*.

## Issues / Contact

- Bug reports: <https://github.com/schlein-lab/variantpaths/issues>
- Schlein Lab homepage: TBD
