#!/bin/bash
# Cross-compile VariantPaths to Windows .exe via cross-rs (Docker).
# Output:  dist/VariantPaths.exe   + dist/igh_3sample.vpf  (demo data)
set -euo pipefail
cd "$(dirname "$0")"

. "$HOME/.cargo/env"

echo "[1/4] regenerating sample .vpf (skip if atlas TSV not present)"
ATLAS_TSV="../phase_d/results/branch_atlas_bubble.tsv"
if [ -f "$ATLAS_TSV" ]; then
    python3 build_vpf.py \
        --tsv "$ATLAS_TSV" \
        --out samples/igh_3sample.vpf \
        --reference GRCh38
else
    echo "  (atlas TSV not found at $ATLAS_TSV — keeping existing samples/igh_3sample.vpf)"
fi

echo "[2/4] cross-compiling x86_64-pc-windows-gnu"
CARGO_TARGET_DIR=target-cross cross build --release --target x86_64-pc-windows-gnu

echo "[3/4] bundling dist/"
mkdir -p dist
cp target-cross/x86_64-pc-windows-gnu/release/variantpaths.exe dist/VariantPaths.exe
cp samples/igh_3sample.vpf dist/
[ -f samples/igh_3sample.vpz ] && cp samples/igh_3sample.vpz dist/

echo "[4/4] done"
echo "  $(du -h dist/VariantPaths.exe | awk '{print $1}')   dist/VariantPaths.exe"
echo "  $(du -h dist/igh_3sample.vpf  | awk '{print $1}')   dist/igh_3sample.vpf"
[ -f dist/igh_3sample.vpz ] && \
    echo "  $(du -h dist/igh_3sample.vpz | awk '{print $1}')   dist/igh_3sample.vpz"
echo
echo "Run on Windows: double-click VariantPaths.exe, or"
echo "  VariantPaths.exe igh_3sample.vpf [igh_3sample.vpz] [reference.fa]"
