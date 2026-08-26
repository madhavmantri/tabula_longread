#!/usr/bin/env python
"""Scatter a CB-sorted subsampled BAM into N shard BAMs in a single pass.

Reads the shard CB lists produced by plan_dedup_shards.py, builds a CB→shard map,
then streams the input BAM once and dispatches each read to its shard's output BAM.
Reads whose CB is not in any shard list (non-real-cells, missing CB tag) are dropped
— matching the default behavior of `isoseq groupdedup` which skips non-real-cells.

Runs `pbindex` on each shard BAM at the end so isoseq groupdedup can consume them.
"""
import argparse
import os
import subprocess
import sys

import pysam


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-bam", required=True)
    p.add_argument("--shard-cb-dir", required=True,
                   help="Directory containing {sid:03d}.cbs.txt files")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--n-shards", type=int, required=True)
    p.add_argument("--cb-tag", default="CB")
    p.add_argument("--threads", type=int, default=4)
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    cb_to_shard = {}
    for sid in range(args.n_shards):
        path = os.path.join(args.shard_cb_dir, f"{sid:03d}.cbs.txt")
        with open(path) as f:
            for line in f:
                cb = line.strip()
                if cb:
                    cb_to_shard[cb] = sid
    print(f"Loaded {len(cb_to_shard)} CBs across {args.n_shards} shards",
          file=sys.stderr)

    infile = pysam.AlignmentFile(args.input_bam, "rb", check_sq=False,
                                 threads=args.threads)
    shard_paths = [os.path.join(args.out_dir, f"{sid:03d}.bam")
                   for sid in range(args.n_shards)]
    outfiles = [pysam.AlignmentFile(p, "wb", template=infile, threads=1)
                for p in shard_paths]
    counts = [0] * args.n_shards
    dropped = 0
    total = 0

    for read in infile.fetch(until_eof=True):
        total += 1
        if not read.has_tag(args.cb_tag):
            dropped += 1
            continue
        sid = cb_to_shard.get(read.get_tag(args.cb_tag))
        if sid is None:
            dropped += 1
            continue
        outfiles[sid].write(read)
        counts[sid] += 1

    for f in outfiles:
        f.close()
    infile.close()

    print(f"Read {total:,} records; wrote {total - dropped:,}; "
          f"dropped {dropped:,} (no CB or non-real-cell)", file=sys.stderr)
    print(f"Per-shard reads: min={min(counts):,} max={max(counts):,}",
          file=sys.stderr)

    for path in shard_paths:
        subprocess.check_call(["pbindex", "-j", str(args.threads), path])


if __name__ == "__main__":
    main()
