#!/usr/bin/env python3
"""
.bbz (VariantPaths sequence Z-archive) builder.

Companion sidecar to a .bbf — contains the actual ALT nucleotide
sequences per bubble.  Loading the .bbz unlocks:
  - sequence display of each ALT in the browser
  - FASTA export per bubble
  - MSA viewer at loci with multiple bubbles

Pipeline:
  1. Read TSV (`branch_atlas_bubble.tsv`) — bubble_name encodes
     entry_node / exit_node like `bubble_3246_n401_416`.
  2. Read the per-sample GFA (S-lines = segment sequences,
     L-lines = links between segments).
  3. For each bubble, find up to N simple paths from entry_node to
     exit_node via DFS (capped at MAX_DEPTH and MAX_PATHS).
  4. Concatenate node sequences along each path.
  5. Write a binary .bbz:

  Header (32 byte LE):
    magic       u32  "VPZ1"
    version     u16  = 1
    flags       u16  bit0 = zstd compressed payload
    n_bubbles   u32
    bbf_built_unix u64    (link to companion .bbf)
    reserved    u32 * 3

  Per-bubble (zstd-compressed payload):
    u32 bubble_name_len
    bytes bubble_name      (UTF-8)
    u8  n_alts
    For each alt:
      u32 path_len_nodes
      For each node: u32 node_id
      u32 alt_seq_len
      bytes alt_seq        (ASCII ACGTN)

Run on the HPC where GFAs live (avoids syncing 100 MB GFAs locally):

  ssh hummel-login
  python3 build_vpz.py \
      --tsv  results/branch_atlas_bubble.tsv \
      --gfa-map "HG00272_LCL=/beegfs/.../HG00272.gfa,HG01530_LCL=/beegfs/.../HG01530.gfa,HG002_Blood=/beegfs/.../HG002_blood.gfa" \
      --max-bubbles-per-sample 200 \
      --max-paths-per-bubble 5 \
      --max-path-depth 600 \
      --max-alt-len 10000 \
      --out  igh_3sample.vpz
"""
from __future__ import annotations
import argparse, csv, re, struct, sys, time
from pathlib import Path
from collections import defaultdict
import zstandard

MAGIC      = b"VPZ1"
VERSION    = 1
FLAG_ZSTD  = 1 << 0
HEADER_SZ  = 32
NODE_RE    = re.compile(r"_n(\d+)_(\d+)$")  # bubble_<idx>_n<entry>_<exit>


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tsv", required=True,
                   help="branch_atlas_bubble.tsv")
    p.add_argument("--gfa-map", required=True,
                   help="comma-separated SAMPLE=PATH list")
    p.add_argument("--out", required=True)
    p.add_argument("--max-bubbles-per-sample", type=int, default=200,
                   help="cap to keep .bbz emailable")
    p.add_argument("--max-paths-per-bubble", type=int, default=5)
    p.add_argument("--max-path-depth", type=int, default=600,
                   help="abort DFS if path grows longer than this many nodes")
    p.add_argument("--max-alt-len", type=int, default=10_000,
                   help="truncate concatenated ALT sequence past this many bp")
    p.add_argument("--no-compress", action="store_true")
    p.add_argument("--bbf-built-unix", type=int, default=0,
                   help="pin the .bbz to a specific .bbf build (optional)")
    return p.parse_args()


def parse_gfa(path: str) -> tuple[dict[int, str], dict[int, list[int]]]:
    """Parse minimal GFA: returns segment_id->sequence, segment_id->[neighbor ids].
    Treats links as undirected (we walk in both directions during DFS)."""
    segs: dict[int, str] = {}
    adj: dict[int, list[int]] = defaultdict(list)
    n_links = 0
    with open(path) as f:
        for line in f:
            if not line:
                continue
            tag = line[0]
            if tag == "S":
                # S\t<id>\t<seq>\t...
                parts = line.split("\t", 3)
                try:
                    sid = int(parts[1])
                except ValueError:
                    continue
                seq = parts[2].rstrip("\n")
                segs[sid] = seq
            elif tag == "L":
                # L\t<from>\t<from_orient>\t<to>\t<to_orient>\t<overlap>...
                parts = line.split("\t", 5)
                if len(parts) < 5:
                    continue
                try:
                    a = int(parts[1]); b = int(parts[3])
                except ValueError:
                    continue
                adj[a].append(b)
                adj[b].append(a)
                n_links += 1
    print(f"  GFA {path}: {len(segs)} segments, {n_links} links", file=sys.stderr)
    return segs, dict(adj)


def find_alt_paths(adj, segs, entry: int, exit_: int,
                   max_paths: int, max_depth: int) -> list[list[int]]:
    """DFS up to max_paths simple paths from entry to exit, no node revisit."""
    if entry not in segs or exit_ not in segs:
        return []
    paths: list[list[int]] = []
    stack: list[tuple[list[int], set[int]]] = [([entry], {entry})]
    while stack and len(paths) < max_paths:
        path, seen = stack.pop()
        if path[-1] == exit_ and len(path) > 1:
            paths.append(path)
            continue
        if len(path) >= max_depth:
            continue
        last = path[-1]
        for nxt in adj.get(last, ()):
            if nxt in seen:
                continue
            stack.append((path + [nxt], seen | {nxt}))
    return paths


def concat_seq(path: list[int], segs: dict[int, str], max_len: int) -> bytes:
    out = bytearray()
    for sid in path:
        s = segs.get(sid, "")
        out.extend(s.encode("ascii"))
        if len(out) >= max_len:
            return bytes(out[:max_len])
    return bytes(out)


def main():
    args = parse_args()

    # ---- load TSV, group bubbles by sample, parse entry/exit nodes ----
    rows = list(csv.DictReader(open(args.tsv), delimiter="\t"))
    by_sample: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        name = r.get("bubble_name", "")
        m = NODE_RE.search(name)
        if not m:
            continue
        r["_entry_node"] = int(m.group(1))
        r["_exit_node"]  = int(m.group(2))
        by_sample[r["sample"]].append(r)

    print(f"  loaded {sum(len(v) for v in by_sample.values())} bubbles from "
          f"{len(by_sample)} samples", file=sys.stderr)

    # ---- per-sample bubble selection: top-K by total_reads then VAF ----
    selected: list[tuple[str, dict, int, int]] = []
    for sample, rs in by_sample.items():
        rs.sort(key=lambda r: (
            -float(r.get("total_reads", 0) or 0),
            float(r.get("min_vaf", 0) or 0),
        ))
        selected.extend((sample, r, r["_entry_node"], r["_exit_node"])
                        for r in rs[:args.max_bubbles_per_sample])
    print(f"  selected {len(selected)} bubbles for sequence extraction "
          f"(cap {args.max_bubbles_per_sample}/sample)", file=sys.stderr)

    # ---- load GFAs ----
    gfas: dict[str, tuple[dict, dict]] = {}
    for kv in args.gfa_map.split(","):
        sample, path = kv.split("=", 1)
        gfas[sample] = parse_gfa(path)

    # ---- extract alt-paths + sequences ----
    output_records: list[tuple[str, list[tuple[list[int], bytes]]]] = []
    n_paths_total = 0
    n_no_path = 0
    for sample, r, entry, exit_ in selected:
        if sample not in gfas:
            print(f"  warning: no GFA for sample {sample}", file=sys.stderr)
            continue
        segs, adj = gfas[sample]
        paths = find_alt_paths(adj, segs, entry, exit_,
                               args.max_paths_per_bubble, args.max_path_depth)
        if not paths:
            n_no_path += 1
            continue
        alts: list[tuple[list[int], bytes]] = []
        for p in paths:
            seq = concat_seq(p, segs, args.max_alt_len)
            alts.append((p, seq))
            n_paths_total += 1
        output_records.append((r["bubble_name"], alts))

    print(f"  extracted {n_paths_total} alt-paths across "
          f"{len(output_records)} bubbles  ({n_no_path} bubbles had no path)",
          file=sys.stderr)

    # ---- serialize body ----
    body = bytearray()
    for bname, alts in output_records:
        nb = bname.encode("utf-8")
        body += struct.pack("<I", len(nb)) + nb
        body += struct.pack("<B", min(255, len(alts)))
        for path, seq in alts:
            body += struct.pack("<I", len(path))
            for sid in path:
                body += struct.pack("<I", sid)
            body += struct.pack("<I", len(seq)) + seq

    flags = 0
    if not args.no_compress:
        body = zstandard.ZstdCompressor(level=19).compress(bytes(body))
        flags |= FLAG_ZSTD

    header = struct.pack("<4sHHIQIII",
        MAGIC, VERSION, flags, len(output_records),
        args.bbf_built_unix or int(time.time()),
        0, 0, 0)
    assert len(header) == HEADER_SZ, len(header)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "wb") as f:
        f.write(header)
        f.write(body)
    sz = Path(args.out).stat().st_size
    print(f"wrote {args.out}  {sz/1024:.1f} KB ({n_paths_total} alt-paths embedded)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
