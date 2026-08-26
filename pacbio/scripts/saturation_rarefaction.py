#!/usr/bin/env python3
"""
Sequencing saturation and rarefaction analysis for PacBio Iso-Seq data.

Reads a pre-dedup sorted BAM, extracts (CB, UMI) pairs, and computes
rarefaction analytically using the formula:
    P(detect molecule at fraction f) = 1 - (1-f)^k
where k is the number of duplicate reads for that molecule.
"""

import argparse
import sys
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute sequencing saturation and rarefaction from a BAM file."
    )
    parser.add_argument("--bam", required=True, help="Input sorted BAM file")
    parser.add_argument("--out-prefix", required=True, help="Output prefix for files")
    parser.add_argument("--barcode-tag", default="CB", help="BAM tag for cell barcode (default: CB)")
    parser.add_argument("--umi-tag", default="XM", help="BAM tag for UMI (default: XM)")
    parser.add_argument("--threads", type=int, default=4, help="Threads for BAM reading (default: 4)")
    return parser.parse_args()


def count_molecules(bam_path, barcode_tag, umi_tag, threads):
    """Read BAM and count reads per unique (barcode, UMI) molecule."""
    molecule_counts = Counter()
    n_reads = 0
    n_skipped = 0

    with pysam.AlignmentFile(bam_path, "rb", check_sq=False, threads=threads) as bam:
        for read in bam:
            if not read.has_tag(barcode_tag) or not read.has_tag(umi_tag):
                n_skipped += 1
                continue
            cb = read.get_tag(barcode_tag)
            umi = read.get_tag(umi_tag)
            if cb == "-" or umi == "-":
                n_skipped += 1
                continue
            molecule_counts[(cb, umi)] += 1
            n_reads += 1

    print(f"Total reads with valid CB+UMI: {n_reads:,}", file=sys.stderr)
    print(f"Reads skipped (missing/invalid tags): {n_skipped:,}", file=sys.stderr)
    print(f"Unique molecules (CB+UMI pairs): {len(molecule_counts):,}", file=sys.stderr)
    return molecule_counts, n_reads


def analytical_rarefaction(molecule_counts, total_reads, fractions):
    """
    Compute expected unique molecules at each subsampling fraction analytically.

    For each molecule with k duplicate reads, the probability of detecting it
    when sampling fraction f of reads is: 1 - (1-f)^k

    Also computes cells detected (unique barcodes with at least one detected molecule).
    """
    # Build a DataFrame of (barcode, k) per molecule — one row per unique molecule
    mol_df = pd.DataFrame(
        [(cb, k) for (cb, _umi), k in molecule_counts.items()],
        columns=["barcode", "k"],
    )

    # --- Bulk rarefaction: group by k only ---
    # Many molecules share the same k, so collapse millions of molecules into ~hundreds of rows
    k_hist = mol_df["k"].value_counts().sort_index()
    k_vals = k_hist.index.values.astype(np.float64)   # unique k values
    k_counts = k_hist.values.astype(np.float64)        # how many molecules at each k

    # --- Per-barcode rarefaction: group by (barcode, k) ---
    bc_k_groups = mol_df.groupby(["barcode", "k"]).size().reset_index(name="n_molecules")
    bc_k_vals = bc_k_groups["k"].values.astype(np.float64)
    bc_k_n = bc_k_groups["n_molecules"].values.astype(np.float64)
    bc_k_barcodes = bc_k_groups["barcode"].values

    results = []
    for frac in fractions:
        read_count = int(round(total_reads * frac))
        q = 1.0 - frac  # probability a single read is NOT sampled

        # Bulk: expected unique molecules = sum(count_k * (1 - q^k))
        detect_by_k = 1.0 - q ** k_vals
        expected_unique = float(np.dot(k_counts, detect_by_k))

        # Saturation
        saturation = 1.0 - (expected_unique / read_count) if read_count > 0 else 0.0

        # Per-barcode metrics using the (barcode, k) groups
        # Expected molecules per (barcode, k) group = n_molecules * (1 - q^k)
        detect_by_bc_k = 1.0 - q ** bc_k_vals
        expected_per_group = bc_k_n * detect_by_bc_k

        # Sum expected molecules per barcode
        bc_expected_df = pd.DataFrame({
            "barcode": bc_k_barcodes,
            "expected_mol": expected_per_group,
        })
        per_bc = bc_expected_df.groupby("barcode")["expected_mol"].sum()

        # Cells detected: barcodes where expected molecules > 0.5
        # (more precisely: sum of P(barcode detected) = 1 - product(q^k_j) per barcode)
        # Since expected_mol = sum(1 - q^k_j) for each barcode, and for detection we need
        # P(at least one molecule) = 1 - product(q^k_j), compute log-sum for accuracy
        log_p_not_detect = bc_k_n * bc_k_vals * np.log(q) if q > 0 else np.full_like(bc_k_vals, -np.inf)
        bc_log_df = pd.DataFrame({
            "barcode": bc_k_barcodes,
            "log_p_not": log_p_not_detect,
        })
        bc_log_sum = bc_log_df.groupby("barcode")["log_p_not"].sum()
        expected_cells = float(np.sum(1.0 - np.exp(bc_log_sum.values)))

        median_mol_per_cell = float(per_bc.median())

        results.append({
            "Percent_reads": int(round(frac * 100)),
            "Read_count": read_count,
            "Unique_molecules": round(expected_unique),
            "Saturation": round(saturation, 6),
            "Cells_detected": round(expected_cells),
            "Median_molecules_per_cell": round(median_mol_per_cell, 1),
        })

    return pd.DataFrame(results)


def per_cell_saturation(molecule_counts):
    """Compute per-cell saturation at 100% of reads."""
    mol_df = pd.DataFrame(
        [(cb, k) for (cb, _umi), k in molecule_counts.items()],
        columns=["Barcode", "k"],
    )
    per_cell = mol_df.groupby("Barcode").agg(
        Total_reads=("k", "sum"),
        Unique_molecules=("k", "size"),
    )
    per_cell["Saturation"] = (1.0 - per_cell["Unique_molecules"] / per_cell["Total_reads"]).round(6)
    return per_cell.sort_values("Total_reads", ascending=False).reset_index()


def plot_saturation(rarefaction_df, per_cell_df, out_pdf):
    """Generate 4-panel saturation/rarefaction figure."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    fig.suptitle("Sequencing Saturation & Rarefaction", fontsize=14, fontweight="bold")

    # (A) Sequencing saturation vs read count
    ax = axes[0, 0]
    ax.plot(rarefaction_df["Read_count"], rarefaction_df["Saturation"] * 100, "o-", color="steelblue")
    ax.set_xlabel("Number of reads")
    ax.set_ylabel("Sequencing saturation (%)")
    ax.set_title("(A) Sequencing Saturation")
    ax.set_ylim(0, 100)
    ax.set_xlim(left=0)

    # (B) Unique molecules vs read count (complexity curve)
    ax = axes[0, 1]
    ax.plot(rarefaction_df["Read_count"], rarefaction_df["Unique_molecules"], "o-", color="darkorange")
    ax.set_xlabel("Number of reads")
    ax.set_ylabel("Unique molecules")
    ax.set_title("(B) Library Complexity")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)

    # (C) Cells detected vs read count
    ax = axes[1, 0]
    ax.plot(rarefaction_df["Read_count"], rarefaction_df["Cells_detected"], "o-", color="seagreen")
    ax.set_xlabel("Number of reads")
    ax.set_ylabel("Cells detected")
    ax.set_title("(C) Cells Detected")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)

    # (D) Per-cell saturation distribution (histogram)
    ax = axes[1, 1]
    ax.hist(per_cell_df["Saturation"] * 100, bins=50, color="mediumpurple", edgecolor="white")
    ax.set_xlabel("Per-cell saturation (%)")
    ax.set_ylabel("Number of cells")
    ax.set_title("(D) Per-cell Saturation Distribution")
    ax.set_xlim(0, 100)

    plt.tight_layout()
    fig.savefig(out_pdf, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved plot: {out_pdf}", file=sys.stderr)


def main():
    args = parse_args()

    print(f"Reading BAM: {args.bam}", file=sys.stderr)
    molecule_counts, total_reads = count_molecules(
        args.bam, args.barcode_tag, args.umi_tag, args.threads
    )

    if total_reads == 0:
        print("ERROR: No valid reads found in BAM.", file=sys.stderr)
        sys.exit(1)

    # Rarefaction at 10% increments
    fractions = [f / 10.0 for f in range(1, 11)]
    print("Computing analytical rarefaction...", file=sys.stderr)
    rarefaction_df = analytical_rarefaction(molecule_counts, total_reads, fractions)

    # Per-cell saturation
    print("Computing per-cell saturation...", file=sys.stderr)
    per_cell_df = per_cell_saturation(molecule_counts)

    # Write outputs
    rarefaction_tsv = f"{args.out_prefix}.rarefaction.tsv"
    per_cell_tsv = f"{args.out_prefix}.per_cell_saturation.tsv"
    pdf_path = f"{args.out_prefix}.saturation_rarefaction.pdf"

    rarefaction_df.to_csv(rarefaction_tsv, sep="\t", index=False)
    print(f"Saved: {rarefaction_tsv}", file=sys.stderr)

    per_cell_df.to_csv(per_cell_tsv, sep="\t", index=False)
    print(f"Saved: {per_cell_tsv}", file=sys.stderr)

    plot_saturation(rarefaction_df, per_cell_df, pdf_path)

    # Print summary
    print("\n--- Rarefaction Summary ---", file=sys.stderr)
    print(rarefaction_df.to_string(index=False), file=sys.stderr)
    print(f"\nTotal cells: {per_cell_df.shape[0]:,}", file=sys.stderr)
    print(f"Median per-cell saturation: {per_cell_df['Saturation'].median():.1%}", file=sys.stderr)


if __name__ == "__main__":
    main()
