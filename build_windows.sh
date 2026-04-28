#!/bin/bash
# Cross-compile VariantPaths to Windows .exe via cross-rs (Docker).
# Output:  dist/VariantPaths.exe   + dist/igh_3sample.vpf  (demo data)
set -euo pipefail
cd "$(dirname "$0")"

. "$HOME/.cargo/env"

echo "[1/4] regenerating sample .bbf"
python3 build_vpf.py \
    --tsv ../results/branch_atlas_bubble.tsv \
    --out samples/igh_3sample.vpf \
    --reference GRCh38

echo "[2/4] cross-compiling x86_64-pc-windows-gnu"
# Isolated target dir avoids host/container glibc-script conflict.
CARGO_TARGET_DIR=target-cross cross build --release --target x86_64-pc-windows-gnu

echo "[3/4] bundling dist/"
mkdir -p dist
cp target-cross/x86_64-pc-windows-gnu/release/variantpaths.exe dist/VariantPaths.exe
cp samples/igh_3sample.vpf dist/

echo "[4/4] done"
echo "  $(du -h dist/VariantPaths.exe | awk '{print $1}')   dist/VariantPaths.exe"
echo "  $(du -h dist/igh_3sample.vpf  | awk '{print $1}')   dist/igh_3sample.vpf"
echo
echo "Run on Windows: double-click VariantPaths.exe, or"
echo "  VariantPaths.exe igh_3sample.vpf"
