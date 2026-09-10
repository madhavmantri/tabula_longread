#!/usr/bin/env python3
"""Differential isoform fraction analysis for PacBio Iso-Seq AnnData.

Identifies isoforms whose usage (fraction) shifts significantly between
cell groups.  Analogous to scanpy's ``rank_genes_groups`` but designed for
bounded compositional isoform-fraction data.

Three statistical methods are supported:

- **wilcoxon** (default): Mann-Whitney U test via ``scipy.stats.mannwhitneyu``.
  Robust for bounded [0, 1] data with many zeros.
- **t-test**: Welch's t-test via ``scipy.stats.ttest_ind``.  Quick exploratory.
- **beta_regression**: Beta regression via ``statsmodels.othermod.betareg.BetaModel``.
  More powerful for proportions; requires statsmodels.

Results are stored in ``adata.uns`` using the same structured-array format as
scanpy's ``rank_genes_groups``.

Usage as a module
-----------------
>>> from differential_isoform_fraction import rank_isoform_fractions
>>> from differential_isoform_fraction import get_rank_isoform_fractions_df
>>> rank_isoform_fractions(adata, groupby="cell_type")
>>> df = get_rank_isoform_fractions_df(adata, group="T_cell")

Usage from CLI
--------------
$ python differential_isoform_fraction.py input.h5ad \\
      --groupby cell_type --method wilcoxon -o output.h5ad
"""

import argparse
import logging
import os
import warnings
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import issparse, csc_matrix
import anndata as ad

logger = logging.getLogger(__name__)

# --- Oversubscription guards (shared by all process-parallel paths) -----------
# Two multiplicative sources of oversubscription must both be bounded:
#   1. the number of worker PROCESSES (n_jobs / n_threads), and
#   2. the BLAS/OpenMP THREADS each worker's numpy/scipy spins up (defaults to
#      every core on the node -> processes x cores threads).
# We clamp (1) to the cores actually granted to this process, and cap (2) to a
# single thread per worker via threadpoolctl. Since the pools use fork(), a
# `threadpool_limits(1)` context around pool creation propagates to the children,
# so total threads = n_workers x 1 <= allocated cores. No kernel/env setup needed.
try:
    from threadpoolctl import threadpool_limits as _threadpool_limits
except ImportError:  # pragma: no cover - fall back to a no-op if unavailable
    from contextlib import contextmanager as _contextmanager

    @_contextmanager
    def _threadpool_limits(limits=None, user_api=None):
        yield


def _available_cpus():
    """CPUs actually usable by this process. Under SLURM the cgroup restricts
    sched_getaffinity() to the allocation, so this never exceeds what the job
    requested (falls back to os.cpu_count() off Linux)."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:  # pragma: no cover - non-Linux
        return max(1, os.cpu_count() or 1)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _benjamini_hochberg(pvals):
    """Apply Benjamini-Hochberg FDR correction.

    Parameters
    ----------
    pvals : numpy.ndarray
        Raw p-values.

    Returns
    -------
    numpy.ndarray
        Adjusted p-values (same length, same order as input).
    """
    n = len(pvals)
    if n == 0:
        return np.array([], dtype=np.float64)

    pvals = np.where(np.isnan(pvals), 1.0, pvals)

    order = np.argsort(pvals)
    pvals_sorted = pvals[order]

    # BH: p_adj[i] = p_sorted[i] * n / rank_i  (rank is 1-based)
    adj = pvals_sorted * n / np.arange(1, n + 1)
    # Enforce monotonicity (from the largest rank downward)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)

    result = np.empty(n, dtype=np.float64)
    result[order] = adj
    return result


def _test_isoform(fracs_group, fracs_ref, method):
    """Run a two-sided test on isoform fractions between two groups.

    Parameters
    ----------
    fracs_group, fracs_ref : numpy.ndarray
        Fraction values for the test and reference cells.
    method : str
        ``'wilcoxon'``, ``'t-test'``, or ``'beta_regression'``.

    Returns
    -------
    tuple[float, float]
        ``(test_statistic, p_value)``.
    """
    if method == "wilcoxon":
        try:
            stat, pval = stats.mannwhitneyu(
                fracs_group, fracs_ref, alternative="two-sided",
            )
        except ValueError:
            stat, pval = 0.0, 1.0
        return float(stat), float(pval)

    if method == "t-test":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            stat, pval = stats.ttest_ind(
                fracs_group, fracs_ref, equal_var=False,
            )
        if np.isnan(pval):
            stat, pval = 0.0, 1.0
        return float(stat), float(pval)

    if method == "beta_regression":
        try:
            from statsmodels.othermod.betareg import BetaModel
        except ImportError:
            raise ImportError(
                "Beta regression requires statsmodels.  "
                "Install with: pip install statsmodels"
            )

        eps = 1e-6
        y = np.concatenate([fracs_group, fracs_ref])
        y = np.clip(y, eps, 1.0 - eps)

        n_g = len(fracs_group)
        n_r = len(fracs_ref)
        x = np.column_stack([
            np.ones(n_g + n_r),
            np.concatenate([np.ones(n_g), np.zeros(n_r)]),
        ])

        try:
            model = BetaModel(y, x)
            result = model.fit(disp=False)
            stat = float(result.tvalues[1])
            pval = float(result.pvalues[1])
        except Exception:
            stat, pval = 0.0, 1.0

        if np.isnan(pval):
            stat, pval = 0.0, 1.0
        return stat, pval

    raise ValueError(f"Unknown method: {method!r}")


def _process_gene_batch_diff(frac_csc, gene_batch, group_to_idx, ref_to_idx,
                              isoform_names, method):
    """Process a batch of genes for differential isoform fraction testing.

    Parameters
    ----------
    frac_csc : scipy.sparse.csc_matrix
        Fraction matrix (cells x isoforms) in CSC format.
    gene_batch : list[tuple[str, numpy.ndarray, list[str]]]
        Each element is ``(gene_name, col_indices, active_groups)``.
    group_to_idx : dict[str, numpy.ndarray]
        Group name -> cell row indices.
    ref_to_idx : dict[str, numpy.ndarray]
        Group name -> reference cell row indices.
    isoform_names : numpy.ndarray
        Isoform names array (length = n_isoforms).
    method : str
        Statistical test method.

    Returns
    -------
    dict[str, dict[str, list]]
        Per-group result accumulators.
    """
    batch_results = {}

    for gene, col_indices_arr, active_groups in gene_batch:
        # Densify only this gene's columns (n_cells x n_isoforms_for_gene)
        sub = frac_csc[:, col_indices_arr].toarray()

        for group in active_groups:
            g_idx = group_to_idx[group]
            r_idx = ref_to_idx[group]

            sub_group = sub[g_idx]
            sub_ref = sub[r_idx]

            # Identify isoforms with any signal
            has_signal = (sub_group.max(axis=0) > 0) | (sub_ref.max(axis=0) > 0)
            if not has_signal.any():
                continue

            signal_local = np.where(has_signal)[0]
            signal_global = col_indices_arr[signal_local]
            sg = sub_group[:, signal_local]
            sr = sub_ref[:, signal_local]
            mean_g = sg.mean(axis=0)
            mean_r = sr.mean(axis=0)
            delta = mean_g - mean_r

            # Vectorized test across all isoforms at once
            if method == "wilcoxon":
                try:
                    stat_arr, pval_arr = stats.mannwhitneyu(
                        sg, sr, alternative="two-sided", axis=0,
                    )
                except ValueError:
                    stat_arr = np.zeros(sg.shape[1])
                    pval_arr = np.ones(sg.shape[1])
                stat_arr = np.atleast_1d(stat_arr).astype(np.float64)
                pval_arr = np.atleast_1d(pval_arr).astype(np.float64)
            elif method == "t-test":
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    stat_arr, pval_arr = stats.ttest_ind(
                        sg, sr, equal_var=False, axis=0,
                    )
                stat_arr = np.atleast_1d(stat_arr).astype(np.float64)
                pval_arr = np.atleast_1d(pval_arr).astype(np.float64)
                nan_mask = np.isnan(pval_arr)
                stat_arr[nan_mask] = 0.0
                pval_arr[nan_mask] = 1.0
            else:
                # beta_regression: must loop per isoform
                n_iso = sg.shape[1]
                stat_arr = np.empty(n_iso, dtype=np.float64)
                pval_arr = np.empty(n_iso, dtype=np.float64)
                for k in range(n_iso):
                    stat_arr[k], pval_arr[k] = _test_isoform(
                        sg[:, k], sr[:, k], method,
                    )

            if group not in batch_results:
                batch_results[group] = {
                    "names": [], "scores": [], "pvals": [],
                    "mean_fraction": [], "mean_fraction_reference": [],
                    "delta_fraction": [], "gene": [],
                }

            acc = batch_results[group]
            acc["names"].extend(isoform_names[signal_global])
            acc["scores"].extend(stat_arr)
            acc["pvals"].extend(pval_arr)
            acc["mean_fraction"].extend(mean_g)
            acc["mean_fraction_reference"].extend(mean_r)
            acc["delta_fraction"].extend(delta)
            acc["gene"].extend([gene] * len(signal_local))

    return batch_results


# ---------------------------------------------------------------------------
# Core API
# ---------------------------------------------------------------------------

def rank_isoform_fractions(
    adata,
    groupby,
    groups="all",
    reference="rest",
    method="wilcoxon",
    layer="isoform_fraction",
    gene_col="gene_names",
    min_cells=10,
    min_gene_frac=0.05,
    corr_method="benjamini-hochberg",
    key_added="rank_isoform_fractions",
    n_threads=4,
    batch_size=500,
):
    """Rank isoforms by differential fraction usage between cell groups.

    For each group, tests whether each isoform's fraction differs
    significantly from the reference.  Results are stored in
    ``adata.uns[key_added]`` in a structured-array format compatible
    with scanpy conventions.

    Parameters
    ----------
    adata : anndata.AnnData
        Must contain ``adata.layers[layer]`` (isoform fractions) and
        ``adata.var[gene_col]`` (isoform-to-gene mapping).
    groupby : str
        Column in ``adata.obs`` defining cell groups.
    groups : ``'all'`` or list[str]
        Which groups to test.
    reference : str
        ``'rest'`` (each group vs all others) or a specific group name.
    method : str
        ``'wilcoxon'`` | ``'t-test'`` | ``'beta_regression'``.
    layer : str
        Layer key for isoform fractions.
    gene_col : str
        ``adata.var`` column mapping isoforms to parent genes.
    min_cells : int
        Minimum cells with the parent gene expressed (nonzero fraction)
        in **both** group and reference to include an isoform.
    min_gene_frac : float
        Minimum fraction of cells expressing the parent gene in **both**
        group and reference.
    corr_method : str
        P-value correction.  Only ``'benjamini-hochberg'`` is supported.
    key_added : str
        Key under which results are stored in ``adata.uns``.
    n_threads : int
        Number of threads for parallel gene-batch processing (default 4).
        Set to 1 to disable threading.  ``scipy.stats.mannwhitneyu`` with
        ``axis=0`` releases the GIL, enabling true parallel speedup.
    batch_size : int
        Number of genes per thread work unit (default 500).

    Returns
    -------
    None
        Modifies ``adata.uns`` in place.
    """
    # --- validation --------------------------------------------------------
    if groupby not in adata.obs.columns:
        raise KeyError(
            f"Column '{groupby}' not found in adata.obs. "
            f"Available: {list(adata.obs.columns)}"
        )
    if layer not in adata.layers:
        raise KeyError(
            f"Layer '{layer}' not found. Available: {list(adata.layers.keys())}. "
            f"Run compute_isoform_fractions() first."
        )
    if gene_col not in adata.var.columns:
        raise KeyError(
            f"Column '{gene_col}' not found in adata.var. "
            f"Available: {list(adata.var.columns)}"
        )
    if corr_method != "benjamini-hochberg":
        raise ValueError(
            f"Only 'benjamini-hochberg' correction is supported, got {corr_method!r}"
        )
    valid_methods = {"wilcoxon", "t-test", "beta_regression"}
    if method not in valid_methods:
        raise ValueError(f"method must be one of {valid_methods}, got {method!r}")

    # --- resolve groups ----------------------------------------------------
    all_groups = adata.obs[groupby].unique()
    if isinstance(groups, str) and groups == "all":
        groups_list = sorted(str(g) for g in all_groups)
    else:
        if isinstance(groups, str):
            groups_list = [groups]
        else:
            groups_list = [str(g) for g in groups]
        missing = set(groups_list) - {str(g) for g in all_groups}
        if missing:
            raise ValueError(f"Groups not found in obs['{groupby}']: {missing}")

    if reference != "rest" and str(reference) not in {str(g) for g in all_groups}:
        raise ValueError(
            f"Reference '{reference}' not found in obs['{groupby}']. "
            f"Available: {sorted(str(g) for g in all_groups)}"
        )

    logger.info(
        "Testing %d groups with method='%s', reference='%s'",
        len(groups_list), method, reference,
    )

    # --- prepare fraction matrix (CSC for column slicing) ------------------
    frac_mat = adata.layers[layer]
    if issparse(frac_mat):
        frac_csc = frac_mat.tocsc()
    else:
        frac_csc = csc_matrix(frac_mat)

    # gene -> list of isoform column indices
    genes = adata.var[gene_col].values
    gene_to_cols = {}
    for idx, g in enumerate(genes):
        gene_to_cols.setdefault(g, []).append(idx)

    isoform_names = np.array(adata.var_names)
    group_labels = adata.obs[groupby].values.astype(str)

    # Pre-compute per-group cell indices
    group_to_idx = {}
    ref_to_idx = {}
    for group in groups_list:
        gmask = group_labels == group
        group_to_idx[group] = np.where(gmask)[0]
        if reference == "rest":
            ref_to_idx[group] = np.where(~gmask)[0]
        else:
            ref_to_idx[group] = np.where(group_labels == str(reference))[0]

    # Pre-filter groups that can never pass min_cells in both sides
    active_groups = [
        g for g in groups_list
        if len(group_to_idx[g]) >= min_cells
        and len(ref_to_idx[g]) >= min_cells
    ]
    n_skipped = len(groups_list) - len(active_groups)
    if n_skipped > 0:
        logger.info(
            "Skipped %d / %d groups with fewer than %d cells in group or reference",
            n_skipped, len(groups_list), min_cells,
        )

    # --- sparse pre-filtering of gene-group pairs ----------------------------
    # Build isoform→gene mapping matrix (n_isoforms x n_unique_genes)
    n_genes_total = len(gene_to_cols)
    gene_names_list = list(gene_to_cols.keys())
    gene_name_to_gidx = {g: i for i, g in enumerate(gene_names_list)}

    map_rows = []
    map_cols = []
    for gene, cols in gene_to_cols.items():
        gidx = gene_name_to_gidx[gene]
        for c in cols:
            map_rows.append(c)
            map_cols.append(gidx)

    isoform_to_gene = csc_matrix(
        (np.ones(len(map_rows), dtype=np.float32),
         (np.array(map_rows, dtype=np.intp),
          np.array(map_cols, dtype=np.intp))),
        shape=(frac_csc.shape[1], n_genes_total),
    )

    # Project binarised fractions to gene space: (n_cells x n_genes)
    frac_binary = frac_csc.copy()
    frac_binary.data = np.ones_like(frac_binary.data)
    gene_expressed = frac_binary @ isoform_to_gene
    gene_expressed.data = np.ones_like(gene_expressed.data)  # binarise
    del frac_binary

    total_expr = np.asarray(gene_expressed.sum(axis=0)).ravel()

    # Precompute reference expression counts when reference is a fixed group
    if reference != "rest" and active_groups:
        ref_idx_common = ref_to_idx[active_groups[0]]
        n_expr_ref_common = np.asarray(
            gene_expressed[ref_idx_common].sum(axis=0)
        ).ravel()
        n_ref_common = len(ref_idx_common)

    active_groups_per_gene = {}  # gene_name -> list[group]
    for group in active_groups:
        g_idx = group_to_idx[group]
        n_group = len(g_idx)
        n_expr_g = np.asarray(gene_expressed[g_idx].sum(axis=0)).ravel()

        if reference == "rest":
            n_expr_r = total_expr - n_expr_g
            n_ref = len(ref_to_idx[group])
        else:
            n_expr_r = n_expr_ref_common
            n_ref = n_ref_common

        passes = (
            (n_expr_g >= min_cells)
            & (n_expr_r >= min_cells)
            & (n_expr_g >= min_gene_frac * n_group)
            & (n_expr_r >= min_gene_frac * n_ref)
        )

        for gi in np.where(passes)[0]:
            gname = gene_names_list[gi]
            if gname not in active_groups_per_gene:
                active_groups_per_gene[gname] = []
            active_groups_per_gene[gname].append(group)

    del gene_expressed, isoform_to_gene

    # Build work items: (gene, col_indices, active_groups_for_gene)
    gene_work_items = []
    for gene, cols in gene_to_cols.items():
        if pd.isna(gene) or str(gene) in ("", "nan"):
            continue
        agroups = active_groups_per_gene.get(gene)
        if not agroups:
            continue
        gene_work_items.append((gene, np.array(cols, dtype=np.intp), agroups))

    logger.info(
        "Pre-filter: %d / %d genes have active group-pairs (skipped %d)",
        len(gene_work_items), n_genes_total, n_genes_total - len(gene_work_items),
    )

    # Initialise per-group result accumulators
    results_per_group = {
        group: {
            "names": [], "scores": [], "pvals": [],
            "mean_fraction": [], "mean_fraction_reference": [],
            "delta_fraction": [], "gene": [],
        }
        for group in groups_list
    }

    # --- batch processing (optionally threaded) ------------------------------
    batches = [
        gene_work_items[i : i + batch_size]
        for i in range(0, len(gene_work_items), batch_size)
    ]
    n_batches = len(batches)

    def _merge_batch_results(batch_results):
        for group, bacc in batch_results.items():
            acc = results_per_group[group]
            for key in acc:
                acc[key].extend(bacc[key])

    if n_threads <= 1 or n_batches <= 1:
        for i, batch in enumerate(batches):
            br = _process_gene_batch_diff(
                frac_csc, batch, group_to_idx, ref_to_idx,
                isoform_names, method,
            )
            _merge_batch_results(br)
            if (i + 1) % 10 == 0 or i == n_batches - 1:
                logger.info("  Processed %d / %d gene batches", i + 1, n_batches)
    else:
        n_threads = max(1, min(n_threads, _available_cpus()))  # never exceed cores
        logger.info("  Using %d threads, %d batches", n_threads, n_batches)
        # One shared process here (threads, not procs), but cap BLAS to 1 so the
        # n_threads workers + BLAS pool stay within the allocation.
        with _threadpool_limits(1), ThreadPoolExecutor(max_workers=n_threads) as pool:
            futures = {
                pool.submit(
                    _process_gene_batch_diff,
                    frac_csc, batch, group_to_idx, ref_to_idx,
                    isoform_names, method,
                ): i
                for i, batch in enumerate(batches)
            }
            done = 0
            for fut in as_completed(futures):
                br = fut.result()
                _merge_batch_results(br)
                done += 1
                if done % 10 == 0 or done == n_batches:
                    logger.info("  Processed %d / %d gene batches", done, n_batches)

    logger.info("  Processed all %d genes with active group-pairs", len(gene_work_items))

    # --- BH correction + sort per group ------------------------------------
    for group in groups_list:
        acc = results_per_group[group]
        pvals_arr = np.array(acc["pvals"], dtype=np.float64)
        pvals_adj = _benjamini_hochberg(pvals_arr)

        order = np.argsort(pvals_adj)
        n_tested = len(order)
        n_sig = int((pvals_adj[order] < 0.05).sum()) if n_tested > 0 else 0
        logger.info(
            "Group '%s': tested %d isoforms, %d significant (adj. p < 0.05)",
            group, n_tested, n_sig,
        )

        results_per_group[group] = {
            "names": np.array(acc["names"], dtype="U100")[order],
            "scores": np.array(acc["scores"], dtype=np.float64)[order],
            "pvals": pvals_arr[order],
            "pvals_adj": pvals_adj[order],
            "mean_fraction": np.array(acc["mean_fraction"], dtype=np.float64)[order],
            "mean_fraction_reference": np.array(acc["mean_fraction_reference"], dtype=np.float64)[order],
            "delta_fraction": np.array(acc["delta_fraction"], dtype=np.float64)[order],
            "gene": np.array(acc["gene"], dtype="U100")[order],
        }

    # --- build structured arrays (scanpy convention) -----------------------
    max_n = max(
        (len(r["names"]) for r in results_per_group.values()), default=0,
    )

    params = {
        "method": method, "groupby": groupby, "reference": reference,
        "corr_method": corr_method, "layer": layer, "gene_col": gene_col,
    }

    if max_n == 0:
        logger.warning("No isoforms passed filters. Check thresholds.")
        adata.uns[key_added] = {"params": params}
        return

    str_dtype = np.dtype([(str(g), "U100") for g in groups_list])
    float_dtype = np.dtype([(str(g), "float64") for g in groups_list])

    uns = {"params": params}

    for key, is_str in [
        ("names", True), ("gene", True),
        ("scores", False), ("pvals", False), ("pvals_adj", False),
        ("mean_fraction", False), ("mean_fraction_reference", False),
        ("delta_fraction", False),
    ]:
        rec = np.zeros(max_n, dtype=str_dtype if is_str else float_dtype)
        for group in groups_list:
            vals = results_per_group[group][key]
            n = len(vals)
            if is_str:
                padded = np.empty(max_n, dtype="U100")
                padded[:] = ""
                if n > 0:
                    padded[:n] = vals
            else:
                padded = np.full(max_n, np.nan)
                if n > 0:
                    padded[:n] = vals
            rec[str(group)] = padded
        uns[key] = rec

    adata.uns[key_added] = uns
    logger.info("Results stored in adata.uns['%s']", key_added)


# ---------------------------------------------------------------------------
# Result extraction
# ---------------------------------------------------------------------------

def _resolve_direction(direction):
    """Normalise a direction argument.

    Accepts ``"both"`` / ``None`` / ``"all"`` (no filter),
    ``"positive"`` / ``"pos"`` / ``"up"`` / ``"+"`` (keep delta > 0),
    ``"negative"`` / ``"neg"`` / ``"down"`` / ``"-"`` (keep delta < 0).
    """
    if direction is None:
        return "both"
    d = str(direction).lower()
    if d in ("both", "all", "two-sided", "twosided"):
        return "both"
    if d in ("positive", "pos", "up", "+"):
        return "positive"
    if d in ("negative", "neg", "down", "-"):
        return "negative"
    raise ValueError(
        f"direction must be one of 'both', 'positive', or 'negative'; got {direction!r}"
    )


def _resolve_sortby(sortby):
    """Normalise a sortby argument for plotting functions.

    Accepted values (case-insensitive):
      ``"pvals_adj"`` / ``"p"`` / ``"pval"`` / ``"pvalue"`` (default ranking),
      ``"delta_fraction"`` / ``"delta"`` / ``"abs_delta_fraction"`` /
      ``"|delta|"`` (rank by |delta_fraction| descending),
      ``"scores"`` / ``"score"`` (rank by |score| descending).
    """
    if sortby is None:
        return "pvals_adj"
    s = str(sortby).lower()
    if s in ("pvals_adj", "pval_adj", "padj", "p", "pval", "pvalue", "pvals"):
        return "pvals_adj"
    if s in ("delta_fraction", "delta", "abs_delta_fraction",
            "abs_delta", "|delta|", "|delta_fraction|"):
        return "delta_fraction"
    if s in ("scores", "score"):
        return "scores"
    raise ValueError(
        f"sortby must be one of 'pvals_adj', 'delta_fraction', or 'scores'; "
        f"got {sortby!r}"
    )


def _reorder_by_sortby(sel_idx, uns, group, sortby):
    """Re-order candidate indices by the chosen sortby metric.

    ``sel_idx`` is the array of candidate positions within
    ``uns[*][group]``.  Returns a new index array in the requested order:

    - ``"pvals_adj"`` -> ascending pvals_adj
    - ``"delta_fraction"`` -> descending ``|delta_fraction|``
    - ``"scores"`` -> descending ``|score|``
    """
    if sel_idx.size == 0:
        return sel_idx
    if sortby == "pvals_adj":
        key = uns["pvals_adj"][group][sel_idx]
        order = np.argsort(key, kind="stable")
    elif sortby == "delta_fraction":
        key = -np.abs(uns["delta_fraction"][group][sel_idx])
        order = np.argsort(key, kind="stable")
    elif sortby == "scores":
        key = -np.abs(uns["scores"][group][sel_idx])
        order = np.argsort(key, kind="stable")
    else:
        raise ValueError(f"Unknown sortby {sortby!r}")
    return sel_idx[order]


def get_rank_isoform_fractions_df(
    adata,
    group=None,
    key="rank_isoform_fractions",
    pval_cutoff=None,
    delta_frac_min=None,
    direction="both",
):
    """Extract differential isoform fraction results as a DataFrame.

    Parameters
    ----------
    adata : anndata.AnnData
        Must contain ``adata.uns[key]`` from a prior call to
        :func:`rank_isoform_fractions`.
    group : str or None
        A specific group name, or ``None`` to return all groups.
    key : str
        ``adata.uns`` key.
    pval_cutoff : float or None
        Filter rows with ``pvals_adj > pval_cutoff``.
    delta_frac_min : float or None
        Filter rows with ``|delta_fraction| < delta_frac_min``.
    direction : str
        Sign filter on ``delta_fraction``.  One of:
        ``"both"`` (default, no filter),
        ``"positive"`` / ``"pos"`` / ``"up"`` (keep ``delta_fraction > 0``),
        ``"negative"`` / ``"neg"`` / ``"down"`` (keep ``delta_fraction < 0``).
        Analogous to scanpy's ``only.pos`` argument.

    Returns
    -------
    pandas.DataFrame
        Columns: group, names, scores, pvals, pvals_adj, mean_fraction,
        mean_fraction_reference, delta_fraction, gene.
    """
    direction = _resolve_direction(direction)
    uns = adata.uns[key]
    field_names = list(uns["names"].dtype.names)

    if group is not None:
        groups_to_show = [str(group)]
        if groups_to_show[0] not in field_names:
            raise KeyError(
                f"Group '{group}' not in results. Available: {field_names}"
            )
    else:
        groups_to_show = field_names

    dfs = []
    for g in groups_to_show:
        names = uns["names"][g]
        mask = names != ""
        if not mask.any():
            continue

        df = pd.DataFrame({
            "group": g,
            "names": names[mask],
            "scores": uns["scores"][g][mask],
            "pvals": uns["pvals"][g][mask],
            "pvals_adj": uns["pvals_adj"][g][mask],
            "mean_fraction": uns["mean_fraction"][g][mask],
            "mean_fraction_reference": uns["mean_fraction_reference"][g][mask],
            "delta_fraction": uns["delta_fraction"][g][mask],
            "gene": uns["gene"][g][mask],
        })

        if pval_cutoff is not None:
            df = df[df["pvals_adj"] <= pval_cutoff]
        if delta_frac_min is not None:
            df = df[df["delta_fraction"].abs() >= delta_frac_min]
        if direction == "positive":
            df = df[df["delta_fraction"] > 0]
        elif direction == "negative":
            df = df[df["delta_fraction"] < 0]

        dfs.append(df)

    if not dfs:
        return pd.DataFrame(
            columns=[
                "group", "names", "scores", "pvals", "pvals_adj",
                "mean_fraction", "mean_fraction_reference",
                "delta_fraction", "gene",
            ]
        )
    return pd.concat(dfs, ignore_index=True)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_rank_isoform_fractions(
    adata,
    groups=None,
    n_isoforms=10,
    key="rank_isoform_fractions",
    ncols=4,
    direction="both",
    sortby="pvals_adj",
):
    """Bar chart of top differential isoforms per group.

    Shows ``delta_fraction`` (group minus reference) for the top
    isoforms selected by the chosen ``sortby`` metric.

    Parameters
    ----------
    adata : anndata.AnnData
        Must contain ``adata.uns[key]``.
    groups : list[str] or None
        Groups to plot.  ``None`` plots all.
    n_isoforms : int
        Number of top isoforms per group.
    key : str
        ``adata.uns`` key.
    ncols : int
        Number of subplot columns.
    direction : str
        Sign filter on ``delta_fraction``.  One of:
        ``"both"`` (default, no filter),
        ``"positive"`` / ``"pos"`` / ``"up"`` (keep ``delta_fraction > 0``),
        ``"negative"`` / ``"neg"`` / ``"down"`` (keep ``delta_fraction < 0``).
    sortby : str
        Ranking metric used to pick the top ``n_isoforms``.  One of:
        ``"pvals_adj"`` (default, ascending adjusted p-value),
        ``"delta_fraction"`` (descending ``|delta_fraction|``),
        ``"scores"`` (descending ``|score|``).

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    direction = _resolve_direction(direction)
    sortby = _resolve_sortby(sortby)

    uns = adata.uns[key]
    field_names = list(uns["names"].dtype.names)

    if groups is not None:
        plot_groups = [str(g) for g in groups]
    else:
        plot_groups = field_names

    nrows = int(np.ceil(len(plot_groups) / ncols))
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize= (2.5 * ncols, 2.0 * nrows),
        squeeze=False,
    )

    for idx, group in enumerate(plot_groups):
        ax = axes[idx // ncols, idx % ncols]
        names = uns["names"][group]
        deltas = uns["delta_fraction"][group]
        mask = names != ""
        if direction == "positive":
            mask = mask & (deltas > 0)
        elif direction == "negative":
            mask = mask & (deltas < 0)

        # Pick the top-N entries passing the direction mask, ordered by sortby.
        cand_idx = np.where(mask)[0]
        cand_idx = _reorder_by_sortby(cand_idx, uns, group, sortby)
        sel_idx = cand_idx[:n_isoforms]
        n_show = sel_idx.size

        if n_show == 0:
            ax.set_title(group)
            ax.text(0.5, 0.5, "no results", transform=ax.transAxes,
                    ha="center", va="center")
            continue

        top_names = names[sel_idx]
        top_delta = deltas[sel_idx]
        top_genes = uns["gene"][group][sel_idx]
        labels = [f"{n}\n({g})" for n, g in zip(top_names, top_genes)]

        colors = ["#d9534f" if d > 0 else "#5bc0de" for d in top_delta]
        ax.barh(np.arange(n_show), top_delta, color=colors)
        ax.set_yticks(np.arange(n_show))
        ax.set_yticklabels(labels, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("delta fraction")
        ax.set_title(group, fontsize=10)
        ax.axvline(0, color="grey", linewidth=0.5, linestyle="--")

    # Hide unused axes
    for idx in range(len(plot_groups), nrows * ncols):
        axes[idx // ncols, idx % ncols].set_visible(False)

    fig.tight_layout()
    return fig


def plot_isoform_fraction_comparison(
    adata,
    group,
    n_isoforms=10,
    key="rank_isoform_fractions",
    direction="both",
    sortby="pvals_adj",
):
    """Grouped bar chart comparing mean fraction in group vs reference.

    Parameters
    ----------
    adata : anndata.AnnData
        Must contain ``adata.uns[key]``.
    group : str
        Group to plot.
    n_isoforms : int
        Number of top isoforms to show.
    key : str
        ``adata.uns`` key.
    direction : str
        Sign filter on ``delta_fraction``.  One of:
        ``"both"`` (default, no filter),
        ``"positive"`` / ``"pos"`` / ``"up"`` (keep ``delta_fraction > 0``),
        ``"negative"`` / ``"neg"`` / ``"down"`` (keep ``delta_fraction < 0``).
    sortby : str
        Ranking metric used to pick the top ``n_isoforms``.  One of:
        ``"pvals_adj"`` (default, ascending adjusted p-value),
        ``"delta_fraction"`` (descending ``|delta_fraction|``),
        ``"scores"`` (descending ``|score|``).

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    direction = _resolve_direction(direction)
    sortby = _resolve_sortby(sortby)

    uns = adata.uns[key]
    group = str(group)
    names = uns["names"][group]
    deltas = uns["delta_fraction"][group]
    mask = names != ""
    if direction == "positive":
        mask = mask & (deltas > 0)
    elif direction == "negative":
        mask = mask & (deltas < 0)

    cand_idx = np.where(mask)[0]
    cand_idx = _reorder_by_sortby(cand_idx, uns, group, sortby)
    sel_idx = cand_idx[:n_isoforms]
    n_show = sel_idx.size

    if n_show == 0:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "no results", transform=ax.transAxes,
                ha="center", va="center")
        return fig

    top_names = names[sel_idx]
    top_genes = uns["gene"][group][sel_idx]
    top_mean = uns["mean_fraction"][group][sel_idx]
    top_mean_ref = uns["mean_fraction_reference"][group][sel_idx]
    top_padj = uns["pvals_adj"][group][sel_idx]

    labels = [f"{n}\n({g})" for n, g in zip(top_names, top_genes)]

    x = np.arange(n_show)
    width = 0.35

    fig, ax = plt.subplots(figsize=(max(5, n_show * 0.5), 3))
    ax.bar(x - width / 2, top_mean, width, label=group, color="#d9534f")
    ref_label = uns["params"].get("reference", "rest")
    ax.bar(x + width / 2, top_mean_ref, width, label=ref_label, color="#5bc0de")

    # Significance stars
    for i, p in enumerate(top_padj):
        if p < 0.001:
            star = "***"
        elif p < 0.01:
            star = "**"
        elif p < 0.05:
            star = "*"
        else:
            star = "ns"
        ymax = max(top_mean[i], top_mean_ref[i])
        ax.text(i, ymax + 0.01, star, ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7, rotation=45, ha="right")
    ax.set_ylabel("Mean isoform fraction")
    ax.set_title(f"Top differential isoforms: {group}")
    ax.legend()
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Gene-level Dirichlet-multinomial isoform-usage test
# ---------------------------------------------------------------------------
#
# WHY: the per-isoform Wilcoxon path (`rank_isoform_fractions`) only tests an
# isoform whose parent gene clears a two-sided min_cells/min_gene_frac floor
# (see lines ~457-462). An isoform that DROPS in the test group sits at low
# fraction there, its gene can fail the group-side floor, and the isoform is
# never tested -> absent from the result -> NaN after a reindex onto the full
# isoform set -> matplotlib silently drops it (non-finite y). The volcano then
# looks one-sided even though `delta_fraction`/`mean_delta` (computed for every
# isoform, tested or not) look balanced. Auditing delta can't reveal this;
# dropping NaN rows can't recover the untested half.
#
# FIX: test each GENE's isoform composition with a Dirichlet-multinomial
# likelihood-ratio test on raw counts -> ONE p-value per gene that jointly
# captures gaining AND losing isoforms, with no per-isoform testability gate.
# The gene p-value is broadcast onto each of the gene's isoforms so a volcano
# (x = per-isoform signed delta_fraction, y = -log10(gene padj)) shows both
# halves. Set include_untested=True to also surface genes that couldn't be
# tested (too few cells) at padj=1 so the coverage gap is visible, not hidden.

def _dm_negloglik_and_grad(theta, C, n_i):
    """Negative Dirichlet-multinomial log-likelihood and gradient wrt theta.

    alpha = exp(theta) (K,). C is (n_cells, K) counts, n_i its row sums.
    logL = sum_i [ gammaln(A) - gammaln(n_i + A)
                   + sum_k (gammaln(C_ik + a_k) - gammaln(a_k)) ], A = sum_k a_k.
    """
    from scipy.special import gammaln, digamma
    alpha = np.exp(theta)
    A = alpha.sum()
    ll = float((gammaln(A) - gammaln(n_i + A)).sum()
               + (gammaln(C + alpha) - gammaln(alpha)).sum())
    g_common = float((digamma(A) - digamma(n_i + A)).sum())
    grad_alpha = g_common + (digamma(C + alpha) - digamma(alpha)).sum(axis=0)
    grad_theta = grad_alpha * alpha          # chain rule: dalpha/dtheta = alpha
    return -ll, -grad_theta


def _dm_fit(C, max_iter=200):
    """MLE of the DM concentration vector for counts C (cells x K).

    Zero-total cells contribute 0 to the likelihood and are dropped.
    Returns (alpha, loglik) or (None, 0.0) on failure / degenerate input.
    """
    from scipy.optimize import minimize
    C = np.ascontiguousarray(C, dtype=np.float64)
    n_i = C.sum(axis=1)
    keep = n_i > 0
    C, n_i = C[keep], n_i[keep]
    K = C.shape[1]
    if C.shape[0] == 0 or K == 0:
        return None, 0.0
    comp = C.sum(axis=0)
    comp = comp / comp.sum() if comp.sum() > 0 else np.full(K, 1.0 / K)
    comp = np.clip(comp, 1e-3, None)
    comp = comp / comp.sum()
    theta0 = np.log(np.clip(comp * K, 1e-6, 1e6))       # precision ~ K
    bounds = [(np.log(1e-8), np.log(1e8))] * K
    try:
        res = minimize(_dm_negloglik_and_grad, theta0, args=(C, n_i),
                       jac=True, method="L-BFGS-B", bounds=bounds,
                       options={"maxiter": max_iter})
        return np.exp(res.x), float(-res.fun)
    except Exception:
        return None, 0.0


def _dm_loglik(alpha, C):
    """DM log-likelihood of counts C under a fixed concentration alpha."""
    from scipy.special import gammaln
    C = np.ascontiguousarray(C, dtype=np.float64)
    n_i = C.sum(axis=1)
    keep = n_i > 0
    C, n_i = C[keep], n_i[keep]
    if C.shape[0] == 0:
        return 0.0
    A = float(alpha.sum())
    return float((gammaln(A) - gammaln(n_i + A)).sum()
                 + (gammaln(C + alpha) - gammaln(alpha)).sum())


def _dm_lrt(Cg, Cr):
    """Dirichlet-multinomial likelihood-ratio test of composition equality.

    H0: group and reference share one concentration vector (K params).
    H1: each side has its own (2K params).  stat = 2*(llH1 - llH0) ~ chi2_K.
    Returns (stat, pval, df).
    """
    K = Cg.shape[1]
    alpha_g, ll_g = _dm_fit(Cg)
    alpha_r, ll_r = _dm_fit(Cr)
    alpha0, _ = _dm_fit(np.vstack([Cg, Cr]))
    if alpha_g is None or alpha_r is None or alpha0 is None:
        return 0.0, 1.0, K
    ll0 = _dm_loglik(alpha0, Cg) + _dm_loglik(alpha0, Cr)
    stat = max(2.0 * ((ll_g + ll_r) - ll0), 0.0)
    return float(stat), float(stats.chi2.sf(stat, K)), K


def _mean_isoform_fraction(C):
    """Mean isoform fraction over cells that express the gene (row sum > 0)."""
    C = np.asarray(C, dtype=np.float64)
    n_i = C.sum(axis=1)
    keep = n_i > 0
    if not keep.any():
        return np.zeros(C.shape[1])
    return (C[keep] / n_i[keep][:, None]).mean(axis=0)


def _cap_isoforms(Cg, Cr, names, max_isoforms):
    """Keep the top (max_isoforms-1) isoforms by total count; lump the rest
    into a single 'other' column for the DM fit (keeps df bounded/stable).
    Returns capped (Cg, Cr) used only for the LRT; per-isoform reporting still
    uses the full signal set."""
    K = Cg.shape[1]
    if K <= max_isoforms:
        return Cg, Cr
    tot = Cg.sum(axis=0) + Cr.sum(axis=0)
    top = np.argsort(tot)[::-1][: max_isoforms - 1]
    rest = np.setdiff1d(np.arange(K), top)
    Cg2 = np.column_stack([Cg[:, top], Cg[:, rest].sum(axis=1)])
    Cr2 = np.column_stack([Cr[:, top], Cr[:, rest].sum(axis=1)])
    return Cg2, Cr2


# Shared read-only state for the DM worker processes. Populated by
# rank_gene_isoform_usage_dm BEFORE the pool is created so that fork()ed workers
# inherit it copy-on-write (no per-task pickling of the count matrix). The MLE
# runs in separate processes -> real parallelism (not GIL-bound like threads).
_DM_STATE = {}


def _dm_gene_worker(item):
    """Test one gene (runs in a worker process). Reads counts from the forked
    `_DM_STATE['Xsub']` (group cells first, then reference cells, split at n_g).
    Returns the same per-gene dict as the serial path (identical numbers)."""
    gene, cols = item
    st = _DM_STATE
    Xsub, n_g = st["Xsub"], st["n_g"]
    dense = np.asarray(Xsub[:, cols].todense())
    Cg_all, Cr_all = dense[:n_g], dense[n_g:]
    sig = (Cg_all.sum(axis=0) + Cr_all.sum(axis=0)) > 0
    if int(sig.sum()) < st["min_isoforms"]:
        return None
    Cg, Cr = Cg_all[:, sig], Cr_all[:, sig]
    names = st["isoform_names"][cols][sig]
    ncg = int((Cg.sum(axis=1) > 0).sum())
    ncr = int((Cr.sum(axis=1) > 0).sum())
    fg, fr = _mean_isoform_fraction(Cg), _mean_isoform_fraction(Cr)
    delta = fg - fr
    tested = (ncg >= st["min_cells"]) and (ncr >= st["min_cells"])
    if tested:
        Cg_fit, Cr_fit = _cap_isoforms(Cg, Cr, names, st["max_isoforms"])
        stat, pval, _ = _dm_lrt(Cg_fit, Cr_fit)
    else:
        if not st["include_untested"]:
            return None
        stat, pval = 0.0, 1.0
    tv = 0.5 * float(np.abs(delta).sum())
    return {
        "gene": gene, "names": names, "n_isoforms": int(sig.sum()),
        "n_cells_group": ncg, "n_cells_ref": ncr,
        "dm_stat": stat, "dm_pval": pval, "tv_distance": tv,
        "mean_fraction": fg, "mean_fraction_reference": fr,
        "delta_fraction": delta, "tested": tested,
    }


def rank_gene_isoform_usage_dm(
    adata,
    groupby,
    group,
    reference="rest",
    counts_layer=None,
    gene_col="gene_name",
    min_cells=10,
    min_isoforms=2,
    max_isoforms=30,
    include_untested=False,
    n_threads=4,
    batch_size=200,
):
    """Gene-level differential isoform usage via a Dirichlet-multinomial LRT.

    Tests one group vs a reference (binary contrast). For each multi-isoform
    gene, fits a Dirichlet-multinomial to the per-cell isoform **counts** and
    runs a likelihood-ratio test of composition equality -> one p-value per
    gene (BH-corrected). Unlike the per-isoform Wilcoxon path, this has no
    per-isoform testability gate, so gaining and losing isoforms of the same
    gene are captured together with no directional bias.

    Parameters
    ----------
    adata : AnnData
        Isoform-level object; ``adata.X`` (or ``counts_layer``) must hold raw
        integer isoform counts, and ``adata.var[gene_col]`` the gene mapping.
    groupby : str
        ``adata.obs`` column with the contrast.
    group : str
        The test group (e.g. the senescent label).
    reference : str
        ``'rest'`` (group vs all other cells) or a specific group label.
    counts_layer : str or None
        Layer with raw counts; ``None`` uses ``adata.X``. NOT the fraction layer.
    gene_col : str
        ``adata.var`` column mapping isoforms to genes.
    min_cells : int
        Minimum cells expressing the gene (row-sum > 0) required in BOTH sides
        for the gene to be tested (a symmetric gene-level floor — it does not
        bias which isoforms are seen within a tested gene).
    min_isoforms : int
        Minimum isoforms (with signal) for a gene to be testable.
    max_isoforms : int
        Cap on isoforms entering the DM fit; extras are lumped into 'other'.
    include_untested : bool
        Also emit rows for multi-isoform genes that failed the cell floor,
        with ``dm_pval=dm_pvals_adj=1.0`` and ``tested=False`` (so a volcano
        shows them at p=1 instead of silently dropping them).
    n_threads : int
        Number of worker **processes** for the per-gene DM fits (name kept for
        back-compat). The MLE is CPU-bound Python, so this uses a fork-based
        ProcessPoolExecutor for real parallelism (threads would be GIL-bound);
        set to 1 for a serial run. Results are identical regardless of value.
    batch_size : int
        Unused by the DM path (kept for signature symmetry with the Wilcoxon API).

    Returns
    -------
    pandas.DataFrame
        One row per (gene isoform) with: group, gene, names, n_isoforms,
        n_cells_group, n_cells_ref, dm_stat, dm_pval, dm_pvals_adj,
        tv_distance (0.5*sum|delta|, gene-level effect size), mean_fraction,
        mean_fraction_reference, delta_fraction, tested. The gene p-value is
        identical across a gene's isoforms; delta_fraction is per-isoform and
        signed, so a volcano (x=delta_fraction, y=-log10(dm_pvals_adj)) shows
        both up- and down-in-group isoforms.
    """
    if groupby not in adata.obs.columns:
        raise KeyError(f"'{groupby}' not in adata.obs")
    if gene_col not in adata.var.columns:
        raise KeyError(f"'{gene_col}' not in adata.var")
    if counts_layer is not None and counts_layer not in adata.layers:
        raise KeyError(f"counts_layer '{counts_layer}' not in adata.layers")

    X = adata.X if counts_layer is None else adata.layers[counts_layer]
    if not issparse(X):
        X = csc_matrix(X)

    genes = adata.var[gene_col].astype(str).values
    isoform_names = np.asarray(adata.var_names)
    labels = adata.obs[groupby].astype(str).values
    g_idx = np.where(labels == str(group))[0]
    if reference == "rest":
        r_idx = np.where(labels != str(group))[0]
    else:
        r_idx = np.where(labels == str(reference))[0]
    if len(g_idx) < min_cells or len(r_idx) < min_cells:
        raise ValueError(
            f"Too few cells: group={len(g_idx)}, reference={len(r_idx)} "
            f"(min_cells={min_cells})"
        )

    # One-time cell slicing: restrict to the group+reference cells ONCE, ordered
    # group-first then reference (split at n_g), so each gene densifies only the
    # relevant rows instead of the whole object. Column-sliceable (CSC) per gene.
    n_g = len(g_idx)
    rows_sel = np.concatenate([g_idx, r_idx])
    Xsub = X.tocsr()[rows_sel].tocsc()

    gene_to_cols = {}
    for i, g in enumerate(genes):
        if g in ("", "nan") or pd.isna(genes[i]):
            continue
        gene_to_cols.setdefault(g, []).append(i)
    work = [
        (g, np.array(cols, dtype=np.intp))
        for g, cols in gene_to_cols.items() if len(cols) >= min_isoforms
    ]
    logger.info(
        "DM isoform-usage: %s vs %s | %d cells vs %d | %d multi-isoform genes",
        group, reference, len(g_idx), len(r_idx), len(work),
    )

    # Publish read-only state, then run each gene's DM fit in a separate process
    # (real parallelism; the per-gene MLE is CPU-bound Python that threads cannot
    # speed up under the GIL). fork() lets workers inherit Xsub copy-on-write, so
    # nothing large is pickled per task. Results are independent of worker order
    # (BH is applied afterward), so this is numerically identical to the serial run.
    _DM_STATE.clear()
    _DM_STATE.update({
        "Xsub": Xsub, "n_g": n_g, "isoform_names": isoform_names,
        "min_isoforms": min_isoforms, "min_cells": min_cells,
        "max_isoforms": max_isoforms, "include_untested": include_untested,
    })
    results = []
    try:
        if n_threads <= 1:
            for item in work:
                r = _dm_gene_worker(item)
                if r is not None:
                    results.append(r)
        else:
            import multiprocessing as _mp
            n_threads = max(1, min(n_threads, _available_cpus()))  # never exceed cores
            chunk = max(1, len(work) // (n_threads * 8) or 1)
            # threadpool_limits(1): forked workers inherit BLAS=1, so total
            # threads = n_threads x 1 <= allocated cores (no oversubscription).
            with _threadpool_limits(1), ProcessPoolExecutor(
                max_workers=n_threads, mp_context=_mp.get_context("fork")
            ) as pool:
                for r in pool.map(_dm_gene_worker, work, chunksize=chunk):
                    if r is not None:
                        results.append(r)
    finally:
        _DM_STATE.clear()

    if not results:
        logger.warning("No genes were testable. Check counts_layer / thresholds.")
        return pd.DataFrame(columns=[
            "group", "gene", "names", "n_isoforms", "n_cells_group",
            "n_cells_ref", "dm_stat", "dm_pval", "dm_pvals_adj", "tv_distance",
            "mean_fraction", "mean_fraction_reference", "delta_fraction", "tested",
        ])

    # BH-correct at the GENE level (one test per gene), then broadcast to isoforms
    tested_mask = np.array([r["tested"] for r in results], dtype=bool)
    gene_pvals = np.array([r["dm_pval"] for r in results], dtype=np.float64)
    gene_padj = np.ones(len(results), dtype=np.float64)
    if tested_mask.any():
        gene_padj[tested_mask] = _benjamini_hochberg(gene_pvals[tested_mask])
    n_sig = int((gene_padj[tested_mask] < 0.05).sum()) if tested_mask.any() else 0
    logger.info(
        "DM isoform-usage: %d genes tested, %d untested (surfaced=%s), "
        "%d significant genes (padj<0.05)",
        int(tested_mask.sum()), int((~tested_mask).sum()), include_untested, n_sig,
    )

    rows = []
    for r, padj in zip(results, gene_padj):
        k = len(r["names"])
        rows.append(pd.DataFrame({
            "group": str(group), "gene": r["gene"], "names": r["names"],
            "n_isoforms": r["n_isoforms"], "n_cells_group": r["n_cells_group"],
            "n_cells_ref": r["n_cells_ref"], "dm_stat": r["dm_stat"],
            "dm_pval": r["dm_pval"], "dm_pvals_adj": padj,
            "tv_distance": r["tv_distance"],
            "mean_fraction": r["mean_fraction"],
            "mean_fraction_reference": r["mean_fraction_reference"],
            "delta_fraction": r["delta_fraction"],
            "tested": np.repeat(r["tested"], k),
        }))
    return pd.concat(rows, ignore_index=True)


def plot_gene_isoform_dm_volcano(
    df, group=None, padj_col="dm_pvals_adj", delta_col="delta_fraction",
    padj_threshold=0.05, delta_threshold=0.05, ax=None, annotate_top=0,
):
    """Per-isoform volcano from :func:`rank_gene_isoform_usage_dm` output.

    x = signed per-isoform ``delta_fraction``, y = -log10(gene ``dm_pvals_adj``).
    Because the y-value is the GENE-level padj broadcast to each isoform, both
    up- and down-in-group isoforms of a significant gene appear (no one-sided
    artefact). Untested isoforms (if surfaced) sit at padj=1 -> y=0.
    """
    import matplotlib.pyplot as plt
    d = df if group is None else df[df["group"] == str(group)]
    d = d.copy()
    padj = np.clip(d[padj_col].to_numpy(dtype=float), 1e-300, 1.0)
    x = d[delta_col].to_numpy(dtype=float)
    y = -np.log10(padj)
    sig = (padj < padj_threshold) & (np.abs(x) >= delta_threshold)
    if ax is None:
        _, ax = plt.subplots(figsize=(4, 4))
    ax.scatter(x[~sig], y[~sig], s=6, c="lightgrey", linewidths=0, label="ns")
    ax.scatter(x[sig & (x > 0)], y[sig & (x > 0)], s=8, c="#d9534f",
               linewidths=0, label="up in group")
    ax.scatter(x[sig & (x < 0)], y[sig & (x < 0)], s=8, c="#5bc0de",
               linewidths=0, label="down in group")
    ax.axhline(-np.log10(padj_threshold), color="grey", lw=0.5, ls="--")
    ax.axvline(delta_threshold, color="grey", lw=0.5, ls="--")
    ax.axvline(-delta_threshold, color="grey", lw=0.5, ls="--")
    ax.set_xlabel("delta isoform fraction (group - reference)")
    ax.set_ylabel("-log10(gene adj. p)")
    if group is not None:
        ax.set_title(str(group))
    ax.legend(fontsize=7, frameon=False)
    if annotate_top and sig.any():
        top = d[sig].reindex(d[sig][padj_col].sort_values().index).head(annotate_top)
        for _, row in top.iterrows():
            ax.annotate(f"{row['names']}\n({row['gene']})",
                        (row[delta_col], -np.log10(max(row[padj_col], 1e-300))),
                        fontsize=6)
    return ax


# ---------------------------------------------------------------------------
# Cell-type-parallel driver (serial per split, splits in parallel)
# ---------------------------------------------------------------------------
#
# rank_gene_isoform_usage_dm's INNER (per-gene) process pool is overhead-bound:
# each gene's MLE is only ~40-65 ms, so fork+IPC costs more than it saves. The
# right granularity is to run each call SERIALLY (fast, no pool) but distribute
# whole cell types (splits) across processes -- few large tasks, negligible
# overhead, near-linear speedup up to the number of splits. Result-identical to
# a serial per-split loop (splits are independent; BH is within each call).

_DM_SPLIT_STATE = {}


def _dm_split_worker(split_val):
    """Run the DM for one split value (one cell type), serially, in a process."""
    st = _DM_SPLIT_STATE
    adata = st["adata"]
    mask = adata.obs[st["split_col"]].astype(str) == str(split_val)
    sub = adata[mask.values]
    lab = sub.obs[st["groupby"]].astype(str)
    n_g = int((lab == str(st["group"])).sum())
    n_r = (int((lab != str(st["group"])).sum()) if st["reference"] == "rest"
           else int((lab == str(st["reference"])).sum()))
    if n_g < st["min_cells"] or n_r < st["min_cells"]:
        return None
    df = rank_gene_isoform_usage_dm(
        sub.copy(), groupby=st["groupby"], group=st["group"],
        reference=st["reference"], counts_layer=st["counts_layer"],
        gene_col=st["gene_col"], min_cells=st["min_cells"],
        min_isoforms=st["min_isoforms"], max_isoforms=st["max_isoforms"],
        include_untested=st["include_untested"], n_threads=1,   # SERIAL inside
    )
    if df is None or not len(df):
        return None
    df[st["split_col"]] = str(split_val)
    return df


def rank_gene_isoform_usage_dm_by_celltype(
    adata, split_col, groupby, group, reference="rest",
    counts_layer=None, gene_col="gene_name", min_cells=10, min_isoforms=2,
    max_isoforms=30, include_untested=False, splits=None, n_jobs=8,
):
    """Gene-level DM run once per value of ``split_col`` (e.g. broad_cell_class),
    each split SERIAL but different splits in PARALLEL processes.

    ``groupby``/``group``/``reference`` must already be set on ``adata`` (e.g.
    create a 'sen_group' column on the full object first). Returns one combined
    DataFrame with a ``split_col`` column; splits with too few cells are skipped.
    Numerically identical to looping the splits serially.
    """
    if groupby not in adata.obs.columns:
        raise KeyError(f"'{groupby}' not in adata.obs")
    if split_col not in adata.obs.columns:
        raise KeyError(f"'{split_col}' not in adata.obs")
    if splits is None:
        splits = sorted(adata.obs[split_col].astype(str).unique())
    else:
        splits = [str(s) for s in splits]

    _DM_SPLIT_STATE.clear()
    _DM_SPLIT_STATE.update({
        "adata": adata, "split_col": split_col, "groupby": groupby,
        "group": group, "reference": reference, "counts_layer": counts_layer,
        "gene_col": gene_col, "min_cells": min_cells, "min_isoforms": min_isoforms,
        "max_isoforms": max_isoforms, "include_untested": include_untested,
    })
    n_jobs = max(1, min(n_jobs, _available_cpus()))  # never spawn more procs than cores
    logger.info("DM by %s: %d splits, %d parallel workers",
                split_col, len(splits), min(n_jobs, len(splits)))
    try:
        if n_jobs <= 1 or len(splits) <= 1:
            out = [_dm_split_worker(s) for s in splits]
        else:
            import multiprocessing as _mp
            # threadpool_limits(1): forked split workers inherit BLAS=1, so total
            # threads = n_workers x 1 <= allocated cores (no oversubscription).
            with _threadpool_limits(1), ProcessPoolExecutor(
                max_workers=min(n_jobs, len(splits)),
                mp_context=_mp.get_context("fork"),
            ) as pool:
                out = list(pool.map(_dm_split_worker, splits))
    finally:
        _DM_SPLIT_STATE.clear()

    out = [d for d in out if d is not None and len(d)]
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Differential isoform fraction analysis."
    )
    parser.add_argument("input", help="Input h5ad file (isoform-level AnnData)")
    parser.add_argument("-o", "--output", required=True, help="Output h5ad file")
    parser.add_argument("--groupby", required=True,
                        help="adata.obs column defining cell groups")
    parser.add_argument("--groups", nargs="+", default="all",
                        help="Groups to test (default: all)")
    parser.add_argument("--reference", default="rest",
                        help="Reference: 'rest' or a group name (default: rest)")
    parser.add_argument("--method", default="wilcoxon",
                        choices=["wilcoxon", "t-test", "beta_regression"],
                        help="Statistical method (default: wilcoxon)")
    parser.add_argument("--layer", default="isoform_fraction",
                        help="Layer key (default: isoform_fraction)")
    parser.add_argument("--gene-col", default="gene_names",
                        help="adata.var column for gene names (default: gene_names)")
    parser.add_argument("--min-cells", type=int, default=10,
                        help="Min cells with gene expressed (default: 10)")
    parser.add_argument("--min-gene-frac", type=float, default=0.05,
                        help="Min fraction of cells with gene expressed (default: 0.05)")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of threads (default: 4)")
    parser.add_argument("--batch-size", type=int, default=500,
                        help="Genes per batch (default: 500)")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
    )

    logger.info("Reading %s", args.input)
    adata = ad.read_h5ad(args.input)
    logger.info("Loaded: %d cells x %d isoforms", *adata.shape)

    rank_isoform_fractions(
        adata,
        groupby=args.groupby,
        groups=args.groups,
        reference=args.reference,
        method=args.method,
        layer=args.layer,
        gene_col=args.gene_col,
        min_cells=args.min_cells,
        min_gene_frac=args.min_gene_frac,
        n_threads=args.threads,
        batch_size=args.batch_size,
    )

    logger.info("Writing %s", args.output)
    adata.write_h5ad(args.output)
    logger.info("Done.")


if __name__ == "__main__":
    main()
