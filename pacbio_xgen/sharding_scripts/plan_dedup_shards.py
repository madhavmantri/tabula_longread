#!/usr/bin/env python
"""Plan CB→shard assignment for parallel isoseq groupdedup.

Reads a tissue's bcstats.tsv, filters to real cells, and bin-packs them across
N shards via longest-processing-time-first (LPT): each cell is assigned to
whichever shard currently has the fewest total reads. This spreads the heavy-tail
barcodes across shards so per-shard wall time is roughly balanced.

Outputs one CB list per shard: {out_dir}/{shard_id:03d}.cbs.txt
"""
import argparse
import csv
import heapq
import os
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bcstats", required=True, help="Path to {tissue}.bcstats.tsv")
    p.add_argument("--n-shards", type=int, required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument(
        "--max-reads-per-barcode",
        type=int,
        default=0,
        help="Cap per-CB workload at this value when bin-packing (matches the "
             "subsample_cell_barcode cap). 0 = uncapped. bcstats reports raw "
             "read counts; if the BAM fed to groupdedup is post-subsample, "
             "planning with raw counts wildly over-weights the deepest CBs.",
    )
    args = p.parse_args()

    cap = args.max_reads_per_barcode
    cells = []
    with open(args.bcstats) as f:
        reader = csv.DictReader(f, delimiter="\t")
        cb_col = reader.fieldnames[0]  # "#BarcodeSequence"
        for row in reader:
            if row.get("RealCell", "").strip() == "cell":
                reads = int(row["NumberOfReads"])
                if cap > 0:
                    reads = min(reads, cap)
                cells.append((reads, row[cb_col].strip()))

    if not cells:
        sys.exit(f"No real cells found in {args.bcstats}")
    if args.n_shards > len(cells):
        sys.exit(f"n_shards={args.n_shards} > number of real cells={len(cells)}")

    cells.sort(reverse=True)  # heaviest first for LPT

    heap = [(0, i) for i in range(args.n_shards)]
    heapq.heapify(heap)
    shard_cbs = [[] for _ in range(args.n_shards)]
    shard_loads = [0] * args.n_shards

    for reads, cb in cells:
        load, sid = heapq.heappop(heap)
        shard_cbs[sid].append(cb)
        shard_loads[sid] = load + reads
        heapq.heappush(heap, (load + reads, sid))

    os.makedirs(args.out_dir, exist_ok=True)
    for sid in range(args.n_shards):
        path = os.path.join(args.out_dir, f"{sid:03d}.cbs.txt")
        with open(path, "w") as f:
            for cb in sorted(shard_cbs[sid]):
                f.write(cb + "\n")

    total = sum(shard_loads)
    print(f"Planned {len(cells)} cells across {args.n_shards} shards", file=sys.stderr)
    print(f"Reads per shard: min={min(shard_loads):,} max={max(shard_loads):,} "
          f"mean={total // args.n_shards:,}", file=sys.stderr)
    print(f"Cells per shard: min={min(len(c) for c in shard_cbs)} "
          f"max={max(len(c) for c in shard_cbs)}", file=sys.stderr)


if __name__ == "__main__":
    main()
