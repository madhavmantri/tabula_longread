#!/usr/bin/env python3
"""Compute per-gene isoform fractions from a PacBio isoform AnnData.

For each gene, the isoform fraction of isoform *i* in cell *c* is defined as:

    IF_{i,c} = count_{i,c} / sum_j(count_{j,c})

where the sum runs over all isoforms *j* belonging to the same gene.
Cells with zero total counts for a gene get fraction 0 for all its isoforms.

The result is stored as a sparse layer in ``adata.layers[layer_name]``.

Usage as a module
-----------------
>>> import scanpy as sc
>>> from isoform_fraction import compute_isoform_fractions
>>> adata = sc.read_h5ad("isoforms.h5ad")
>>> compute_isoform_fractions(adata, gene_col="gene_name")
>>> adata.layers["isoform_fraction"]   # sparse CSR matrix

Usage from CLI
--------------
$ python isoform_fraction.py input.h5ad -o output.h5ad --gene-col gene_name --threads 8
"""

import argparse
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
from scipy.sparse import csr_matrix, lil_matrix
import anndata as ad

logger = logging.getLogger(__name__)


def _process_gene_batch(X_csc, gene_groups_batch):
    """Compute isoform fractions for a batch of genes.

    Parameters
    ----------
    X_csc : scipy.sparse.csc_matrix
        Full expression matrix in CSC format (cells x isoforms).
    gene_groups_batch : list[numpy.ndarray]
        Each element is an array of column indices for one gene.

    Returns
    -------
    list[tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]]
        For each gene: (col_indices, row_indices_of_nonzero, fraction_values).
    """
    results = []
    for col_indices in gene_groups_batch:
        if len(col_indices) == 1:
            # Single-isoform gene: fraction is 1 wherever count > 0.
            col = X_csc[:, col_indices[0]]
            nz_rows = col.nonzero()[0]
            vals = np.ones(len(nz_rows), dtype=np.float32)
            results.append((col_indices, nz_rows, vals))
            continue

        sub = X_csc[:, col_indices]  # cells x n_isoforms_in_gene, CSC
        gene_totals = np.asarray(sub.sum(axis=1)).ravel()  # (n_cells,)

        sub_csr = sub.tocsr()
        nz_rows, nz_local_cols = sub_csr.nonzero()
        nz_vals = np.asarray(sub_csr[nz_rows, nz_local_cols]).ravel()

        row_totals = gene_totals[nz_rows]
        fractions = np.divide(
            nz_vals, row_totals, out=np.zeros_like(nz_vals, dtype=np.float32),
            where=row_totals != 0,
        ).astype(np.float32)

        results.append((col_indices, nz_rows, nz_local_cols, fractions))

    return results


def compute_isoform_fractions(
    adata,
    gene_col="gene_name",
    layer_name="isoform_fraction",
    n_threads=1,
    batch_size=500,
):
    """Compute per-gene isoform fractions and store as a sparse layer.

    Parameters
    ----------
    adata : anndata.AnnData
        Isoform-level AnnData. ``adata.X`` should contain raw counts
        (sparse or dense). ``adata.var[gene_col]`` maps isoforms to genes.
    gene_col : str
        Column in ``adata.var`` that holds the gene name / ID.
    layer_name : str
        Key under which to store the result in ``adata.layers``.
    n_threads : int
        Number of threads for parallel gene-batch processing.
    batch_size : int
        Number of genes per thread work unit.

    Returns
    -------
    None
        Modifies ``adata`` in place by adding ``adata.layers[layer_name]``.
    """
    if gene_col not in adata.var.columns:
        raise KeyError(f"Column '{gene_col}' not found in adata.var. "
                       f"Available: {list(adata.var.columns)}")

    n_cells, n_isoforms = adata.shape
    logger.info("Computing isoform fractions: %d cells x %d isoforms", n_cells, n_isoforms)

    # CSC is efficient for column-slicing (selecting isoforms per gene).
    from scipy.sparse import issparse
    if issparse(adata.X):
        X_csc = adata.X.tocsc()
    else:
        X_csc = csr_matrix(adata.X).tocsc()

    # Group isoform column-indices by gene.
    genes = adata.var[gene_col].values
    gene_to_cols = {}
    for idx, g in enumerate(genes):
        gene_to_cols.setdefault(g, []).append(idx)

    gene_groups = [np.array(cols, dtype=np.intp) for cols in gene_to_cols.values()]
    n_genes = len(gene_groups)
    logger.info("Grouped into %d genes (threads=%d, batch=%d)", n_genes, n_threads, batch_size)

    # Split genes into batches for parallel processing.
    batches = [gene_groups[i:i + batch_size] for i in range(0, n_genes, batch_size)]

    # Allocate output in LIL format (efficient for incremental row/col assignment).
    frac = lil_matrix((n_cells, n_isoforms), dtype=np.float32)

    def _fill_results(batch_results):
        for item in batch_results:
            if len(item) == 3:
                # Single-isoform gene.
                col_indices, nz_rows, vals = item
                col = col_indices[0]
                for r, v in zip(nz_rows, vals):
                    frac[r, col] = v
            else:
                col_indices, nz_rows, nz_local_cols, fractions = item
                for r, lc, v in zip(nz_rows, nz_local_cols, fractions):
                    frac[r, col_indices[lc]] = v

    if n_threads <= 1:
        for i, batch in enumerate(batches):
            results = _process_gene_batch(X_csc, batch)
            _fill_results(results)
            if (i + 1) % 20 == 0 or i == len(batches) - 1:
                logger.info("Processed %d / %d gene batches", i + 1, len(batches))
    else:
        with ThreadPoolExecutor(max_workers=n_threads) as pool:
            futures = {pool.submit(_process_gene_batch, X_csc, batch): i
                       for i, batch in enumerate(batches)}
            done = 0
            for fut in as_completed(futures):
                results = fut.result()
                _fill_results(results)
                done += 1
                if done % 20 == 0 or done == len(batches):
                    logger.info("Processed %d / %d gene batches", done, len(batches))

    adata.layers[layer_name] = frac.tocsr()
    logger.info("Stored layer '%s' (%d nnz)", layer_name, adata.layers[layer_name].nnz)


def main():
    parser = argparse.ArgumentParser(
        description="Compute per-gene isoform fractions and save as an AnnData layer."
    )
    parser.add_argument("input", help="Input h5ad file (isoform-level AnnData)")
    parser.add_argument("-o", "--output", required=True, help="Output h5ad file")
    parser.add_argument("--gene-col", default="gene_name",
                        help="adata.var column mapping isoforms to genes (default: gene_name)")
    parser.add_argument("--layer-name", default="isoform_fraction",
                        help="Layer name for the result (default: isoform_fraction)")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of threads (default: 1)")
    parser.add_argument("--batch-size", type=int, default=500,
                        help="Genes per batch (default: 500)")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    logger.info("Reading %s", args.input)
    adata = ad.read_h5ad(args.input)
    logger.info("Loaded: %d cells x %d isoforms", *adata.shape)

    compute_isoform_fractions(
        adata,
        gene_col=args.gene_col,
        layer_name=args.layer_name,
        n_threads=args.threads,
        batch_size=args.batch_size,
    )

    logger.info("Writing %s", args.output)
    adata.write_h5ad(args.output)
    logger.info("Done.")


if __name__ == "__main__":
    main()
