#!/usr/bin/env python3
"""
TSV -> .bbf (VariantPaths File) converter.

Reads either:
  - branch_atlas_bubble.tsv (per-bubble; preferred)
  - any per-sample bubbles.bed augmented with a 'sample' column

Writes a compact binary file:

  Header (64 byte LE):
    magic        u32  b"VPF1"
    version      u16  = 1
    flags        u16  bit0=zstd_payload
    n_chroms     u32
    n_samples    u32
    n_classes    u32
    n_bubbles    u64
    built_unix   u64
    reference_id [16]byte ASCII NUL-padded
    reserved     u8[12]

  Body (zstd-compressed if flags bit0 set):
    sample_pool, class_pool, chrom_pool, bubble_name_pool, dbvar_id_pool
      each: u32 size_bytes; u32 n_entries; [u32 offset_into_data]*n; bytes
    chrom_index: per chrom: u32 chrom_name_idx, u32 chrom_length_bp,
                            u64 bubble_offset, u32 bubble_count
    bubble_array: 32 bytes/record (little-endian)
      i32 start_bp        bytes 0-3
      i32 end_bp          bytes 4-7
      u16 vaf_q16         bytes 8-9     (vaf*65535)
      u16 chrom_idx       bytes 10-11
      u8  length_log10_q  byte  12      (log10(length+1)*16)
      u8  flags           byte  13
      u8  sample_idx      byte  14
      u8  class_idx       byte  15
      u8  n_alts          byte  16
      u8  total_reads_log byte  17      (log2(R)*8 quantized)
      u8  dbvar_recip_q8  byte  18      (recip*255)
      u8  _pad            byte  19
      u32 bubble_name_idx bytes 20-23   (0xFFFFFFFF = none)
      u32 dbvar_id_idx    bytes 24-27   (0xFFFFFFFF = none)
      u32 reserved        bytes 28-31

Sortierung: bubble records sortiert nach (chrom_idx, start_bp).
"""
from __future__ import annotations
import argparse, csv, io, math, struct, sys, time, zstandard
from pathlib import Path
from collections import defaultdict

MAGIC      = b"VPF1"
VERSION    = 1
FLAG_ZSTD  = 1 << 0
RECORD_SZ  = 32
HEADER_SZ  = 64
NIL_IDX    = 0xFFFFFFFF


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tsv", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--reference", default="GRCh38",
                   help="reference assembly id, freitext (max 16 ASCII bytes)")
    p.add_argument("--no-compress", action="store_true",
                   help="skip zstd compression of the body")
    p.add_argument("--inspect", action="store_true",
                   help="after writing, also read back and report stats")
    return p.parse_args()


def build_pool(strings):
    """Pool layout:
        u32 size_bytes_of_data    (sum of len(s) + 1 for NUL)
        u32 n_entries
        [u32 offset]*n_entries     (offset into data block)
        data bytes (each NUL-terminated UTF-8)
    Returns (pool_bytes, idx_map)."""
    idx_map = {}
    data = bytearray()
    offsets = []
    for s in strings:
        if s in idx_map:
            continue
        idx_map[s] = len(offsets)
        offsets.append(len(data))
        data.extend(s.encode("utf-8"))
        data.append(0)
    buf = bytearray()
    buf += struct.pack("<II", len(data), len(offsets))
    for off in offsets:
        buf += struct.pack("<I", off)
    buf += data
    return bytes(buf), idx_map


def vaf_to_q16(v: float) -> int:
    if v < 0.0: v = 0.0
    if v > 1.0: v = 1.0
    return int(round(v * 65535))


def length_log10_q8(L: int) -> int:
    # log10(L+1) * 16 quantized to u8, max 255 ~ 1e16 (lots of headroom)
    if L < 0: L = 0
    q = int(round(math.log10(L + 1) * 16))
    return min(255, max(0, q))


def total_reads_log_q8(R: int) -> int:
    if R < 1: return 0
    q = int(round(math.log2(R) * 8))   # log2 * 8, max 255 ~ 2^32
    return min(255, max(0, q))


def main():
    args = parse_args()

    rows = list(csv.DictReader(open(args.tsv), delimiter="\t"))
    if not rows:
        sys.exit("empty TSV")
    print(f"  {len(rows)} TSV rows", file=sys.stderr)

    samples, classes, chroms = [], [], []
    sample_idx, class_idx, chrom_idx = {}, {}, {}
    bubble_names, dbvar_ids = [], []

    parsed = []
    n_bad = 0
    for r in rows:
        try:
            entry = int(r["entry_grch38"])
            exit_ = int(r["exit_grch38"])
        except (KeyError, ValueError):
            n_bad += 1; continue
        if entry > exit_:
            entry, exit_ = exit_, entry
        try:    vaf   = float(r.get("min_vaf", "0") or 0)
        except: vaf   = 0.0
        try:    length = int(r.get("length_bp", str(exit_ - entry)) or (exit_ - entry))
        except: length = exit_ - entry
        try:    n_alts = int(r.get("n_alts", "0") or 0)
        except: n_alts = 0
        try:    total_reads = int(r.get("total_reads", "0") or 0)
        except: total_reads = 0
        try:    dbvar_recip = float(r.get("dbvar_recip", "0") or 0)
        except: dbvar_recip = 0.0
        is_shared = (r.get("is_shared", "0") or "0").strip() in ("1","true","TRUE","yes")
        sample = r.get("sample", "unknown")
        cls    = r.get("classification", "UNCLASSIFIED")
        chrom  = r.get("chrom", "chr14")  # IGH default
        bname  = r.get("bubble_name", "")
        dbid   = r.get("dbvar_top_id", "") or r.get("dbvar_id", "")

        if sample not in sample_idx:
            sample_idx[sample] = len(samples); samples.append(sample)
        if cls not in class_idx:
            class_idx[cls] = len(classes); classes.append(cls)
        if chrom not in chrom_idx:
            chrom_idx[chrom] = len(chroms); chroms.append(chrom)

        parsed.append({
            "chrom_idx": chrom_idx[chrom],
            "start": entry, "end": exit_, "vaf": vaf,
            "length": length, "sample_idx": sample_idx[sample],
            "class_idx": class_idx[cls], "n_alts": n_alts,
            "total_reads": total_reads, "dbvar_recip": dbvar_recip,
            "name": bname, "dbvar_id": dbid,
            "is_shared": is_shared,
        })
    print(f"  {n_bad} skipped, kept {len(parsed)}", file=sys.stderr)

    # Sort by (chrom_idx, start)
    parsed.sort(key=lambda b: (b["chrom_idx"], b["start"]))

    # Per-chrom bubble_offset / count + chrom_length_bp
    chrom_count   = [0] * len(chroms)
    chrom_offset  = [0] * len(chroms)
    chrom_length  = [0] * len(chroms)
    last_chrom = -1
    for i, b in enumerate(parsed):
        c = b["chrom_idx"]
        if c != last_chrom:
            chrom_offset[c] = i
            last_chrom = c
        chrom_count[c] += 1
        if b["end"] > chrom_length[c]:
            chrom_length[c] = b["end"]

    # Build pools
    sample_pool, _ = build_pool(samples)
    class_pool,  _ = build_pool(classes)
    chrom_pool,  _ = build_pool(chroms)
    bubble_pool, bn_map = build_pool([b["name"] for b in parsed if b["name"]])
    dbvar_pool,  dv_map = build_pool([b["dbvar_id"] for b in parsed if b["dbvar_id"]])

    # chrom index
    chrom_index_bytes = bytearray()
    for ci in range(len(chroms)):
        chrom_index_bytes += struct.pack("<IIQI",
            ci, chrom_length[ci], chrom_offset[ci], chrom_count[ci])

    # bubble records
    rec_bytes = bytearray()
    fmt = "<iiHHBBBBBBBBIII"
    rec_struct = struct.Struct(fmt)
    assert rec_struct.size == RECORD_SZ, f"record fmt size {rec_struct.size} != {RECORD_SZ}"
    for b in parsed:
        bn_idx = bn_map.get(b["name"], NIL_IDX) if b["name"] else NIL_IDX
        dv_idx = dv_map.get(b["dbvar_id"], NIL_IDX) if b["dbvar_id"] else NIL_IDX
        bubble_flags = 0
        if b.get("is_shared"): bubble_flags |= 0x01  # bit 0 = is_shared
        rec_bytes += rec_struct.pack(
            b["start"], b["end"],
            vaf_to_q16(b["vaf"]),
            b["chrom_idx"] & 0xFFFF,
            length_log10_q8(b["length"]),
            bubble_flags,
            b["sample_idx"] & 0xFF,
            b["class_idx"]  & 0xFF,
            b["n_alts"] & 0xFF,
            total_reads_log_q8(b["total_reads"]),
            int(round(min(1.0, max(0.0, b["dbvar_recip"])) * 255)),
            0,  # _pad
            bn_idx, dv_idx,
            0,  # reserved u32
        )
    assert len(rec_bytes) == len(parsed) * RECORD_SZ, \
        f"record size mismatch: {len(rec_bytes)} vs {len(parsed) * RECORD_SZ}"

    # Body assembly
    body = (sample_pool + class_pool + chrom_pool + bubble_pool + dbvar_pool +
            bytes(chrom_index_bytes) + bytes(rec_bytes))
    flags = 0
    if not args.no_compress:
        cctx = zstandard.ZstdCompressor(level=19)
        body = cctx.compress(body)
        flags |= FLAG_ZSTD

    # Header
    ref_bytes = args.reference.encode("ascii")[:16].ljust(16, b"\x00")
    header = struct.pack("<4sHHIIIQQ16s12s",
        MAGIC, VERSION, flags,
        len(chroms), len(samples), len(classes),
        len(parsed),
        int(time.time()),
        ref_bytes,
        b"\x00" * 12,
    )
    assert len(header) == HEADER_SZ, len(header)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "wb") as f:
        f.write(header)
        f.write(body)

    sz = Path(args.out).stat().st_size
    raw_sz = HEADER_SZ + (len(sample_pool) + len(class_pool) + len(chrom_pool)
                         + len(bubble_pool) + len(dbvar_pool)
                         + len(chrom_index_bytes) + len(rec_bytes))
    print(f"wrote {args.out}  {sz/1024:.1f} KB  ({raw_sz/1024:.1f} KB raw, "
          f"{100*sz/raw_sz:.1f}% size)", file=sys.stderr)
    print(f"  n_bubbles={len(parsed)}  chroms={len(chroms)}  "
          f"samples={len(samples)}  classes={len(classes)}", file=sys.stderr)
    print(f"  reference_id={args.reference!r}", file=sys.stderr)

    if args.inspect:
        # quick read-back to sanity-check
        with open(args.out, "rb") as f:
            data = f.read()
        head = struct.unpack_from("<4sHHIIIQQ16s12s", data, 0)
        print(f"  read back: magic={head[0]!r} version={head[1]} flags={head[2]} "
              f"n_bubbles={head[6]}  ref={head[8].rstrip(bytes([0])).decode()}",
              file=sys.stderr)

if __name__ == "__main__":
    main()
