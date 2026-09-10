#!/usr/bin/env python3
"""Add pigeon/SQANTI classification metadata to an isoform-level AnnData.

Reads a pigeon ``collapse_classification.filtered_lite_classification.txt``
file and merges selected columns into ``adata.var``.

The join is performed on the PacBio isoform ID (e.g. ``PB.3.5``).  When the
AnnData var index uses the pigeon Seurat format ``PB.3.5:GENE``, the gene
suffix is automatically stripped before matching.

Usage as a module
-----------------
>>> import scanpy as sc
>>> from add_classification import add_classification
>>> adata = sc.read_h5ad("isoforms.h5ad")
>>> add_classification(adata, "collapse_classification.filtered_lite_classification.txt")
>>> adata.var["structural_category"]

Usage from CLI
--------------
$ python add_classification.py input.h5ad classification.txt -o output.h5ad
"""

import argparse
import logging

import pandas as pd
import anndata as ad

logger = logging.getLogger(__name__)

DEFAULT_COLUMNS = [
    "structural_category",
    "associated_gene",
    "associated_transcript",
    "subcategory",
    "chrom",
    "strand",
    "length",
    "exons",
    "coding",
    "FSM_class",
    "predicted_NMD",
    "FL",
    "ref_length",
    "ref_exons",
    # SQANTI3 RulesFilter call (Isoform vs Artifact). Only present in the
    # RulesFilter-derived SQANTI3 CSVs; absent from pigeon classification and
    # the QC-derived files, where it is skipped with a warning (see below).
    "filter_result",
]


def _extract_pb_id(var_index):
    """Extract the PB isoform ID from the var index.

    Handles both plain IDs (``PB.3.5``) and Seurat-style IDs
    (``PB.3.5:GENE``).
    """
    return var_index.str.split(":").str[0]


def add_classification(
    adata,
    classification_path,
    columns=None,
    isoform_col=None,
    prefix=None,
    sep="\t",
):
    """Merge pigeon classification columns into ``adata.var``.

    Parameters
    ----------
    adata : anndata.AnnData
        Isoform-level AnnData whose var index (or *isoform_col*) contains
        PacBio isoform IDs.
    classification_path : str or pathlib.Path
        Path to the pigeon classification TSV
        (``collapse_classification.filtered_lite_classification.txt``)
        or a CSV (e.g. combined SQANTI3 classification).
    columns : list[str] or None
        Classification columns to add.  ``None`` uses a sensible default set.
    isoform_col : str or None
        Column in ``adata.var`` that holds the isoform ID.  When ``None``,
        the var *index* is used.
    prefix : str or None
        If provided, prepend this string to each added column name
        (e.g. ``prefix="pigeon_"`` produces ``pigeon_structural_category``).
    sep : str
        Column separator for the classification file.  Default ``"\\t"``
        (TSV).  Use ``","`` for CSV files.

    Returns
    -------
    None
        Modifies ``adata.var`` in place.
    """
    if columns is None:
        columns = list(DEFAULT_COLUMNS)

    # --- read classification ---
    logger.info("Reading classification from %s", classification_path)
    cls_df = pd.read_csv(
        classification_path,
        sep=sep,
        dtype=str,
        usecols=lambda c: c == "isoform" or c in columns,
    )
    logger.info("Classification table: %d isoforms, columns %s",
                len(cls_df), list(cls_df.columns))

    # Keep only requested columns that actually exist in the file.
    available = [c for c in columns if c in cls_df.columns]
    missing = set(columns) - set(available)
    if missing:
        logger.warning("Columns not found in classification file: %s", missing)
    cls_df = cls_df[["isoform"] + available]
    cls_df = cls_df.drop_duplicates(subset="isoform")

    # --- build join key from adata ---
    if isoform_col is not None:
        if isoform_col not in adata.var.columns:
            raise KeyError(
                f"Column '{isoform_col}' not found in adata.var. "
                f"Available: {list(adata.var.columns)}"
            )
        pb_ids = adata.var[isoform_col].astype(str)
    else:
        pb_ids = adata.var.index.to_series().astype(str)

    # Strip ":GENE" suffix if present (pigeon Seurat format).
    if pb_ids.str.contains(":").any():
        pb_ids = _extract_pb_id(pb_ids)

    adata.var["_pb_id"] = pb_ids.values

    # --- merge ---
    cls_df = cls_df.set_index("isoform")
    n_before = adata.var["_pb_id"].isin(cls_df.index).sum()
    logger.info("Matched %d / %d isoforms to classification", n_before, adata.n_vars)

    merged = adata.var[["_pb_id"]].join(cls_df, on="_pb_id", how="left")

    for col in available:
        dest_col = f"{prefix}{col}" if prefix else col
        adata.var[dest_col] = merged[col].values

    adata.var.drop(columns=["_pb_id"], inplace=True)

    first_dest = f"{prefix}{available[0]}" if prefix else available[0]
    n_annotated = adata.var[first_dest].notna().sum() if available else 0
    logger.info("Added %d columns; %d / %d isoforms annotated",
                len(available), n_annotated, adata.n_vars)


def main():
    parser = argparse.ArgumentParser(
        description="Add pigeon classification metadata to an isoform AnnData."
    )
    parser.add_argument("input", help="Input h5ad file (isoform-level AnnData)")
    parser.add_argument("classification",
                        help="Pigeon classification TSV file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output h5ad file")
    parser.add_argument("--columns", nargs="+", default=None,
                        help="Classification columns to add (default: sensible set)")
    parser.add_argument("--isoform-col", default=None,
                        help="adata.var column with isoform IDs (default: var index)")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")

    logger.info("Reading %s", args.input)
    adata = ad.read_h5ad(args.input)
    logger.info("Loaded: %d cells x %d isoforms", *adata.shape)

    add_classification(
        adata,
        args.classification,
        columns=args.columns,
        isoform_col=args.isoform_col,
    )

    logger.info("Writing %s", args.output)
    adata.write_h5ad(args.output)
    logger.info("Done.")


if __name__ == "__main__":
    main()
