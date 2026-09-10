from __future__ import annotations

import gc
import logging
import os
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Literal, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from joblib import Parallel, delayed
from scipy import sparse as sp_sparse
from scipy import stats as scipy_stats
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)

GeneList = List[str]
GeneZScores = List[Tuple[str, float]]
GeneLFC = List[Tuple[str, float]]
MarkerResult = Optional[Union[GeneList, GeneZScores, GeneLFC, pd.DataFrame, dict]]

_OLS_GENE_CHUNK = 2000


def _resolve_n_jobs(n_jobs: int) -> int:
    if n_jobs < 0:
        return max(1, os.cpu_count() + 1 + n_jobs)
    return max(1, n_jobs)


def filter_markers(
    markers: pd.DataFrame,
    *,
    only_positive: bool = True,
    pval_thresh: float = 0.01,
    min_log_fc: float = 0.25,
    min_pct_nz_group: float = 0.0,
    min_pct_nz_difference: float = 0.0,
) -> pd.DataFrame:
    pval_mask = markers["pvals_adj"] < pval_thresh
    pct_mask = markers["pct_nz_group"] > min_pct_nz_group

    if only_positive:
        fc_mask = markers["logfoldchanges"] > min_log_fc
    else:
        fc_mask = markers["logfoldchanges"].abs() > min_log_fc

    markers = markers.loc[pval_mask & fc_mask & pct_mask].copy()
    markers.sort_values("logfoldchanges", ascending=False, inplace=True)

    markers["pct_nz_difference"] = (
        markers["pct_nz_group"] - markers["pct_nz_reference"]
    )

    if min_pct_nz_difference > 0.0:
        drop_mask = (markers["pct_nz_reference"] < 0.01) & (
            markers["pct_nz_difference"] < min_pct_nz_difference
        )
        markers = markers.loc[~drop_mask]

    markers["z"] = markers["logfoldchanges"] * (
        -np.log10(markers["pvals_adj"] + 1e-9)
    )

    return markers


def filter_lm_results(
    de_results: pd.DataFrame,
    *,
    only_positive: bool = True,
    pval_thresh: float = 0.01,
    min_log_fc: float = 0.25,
) -> pd.DataFrame:
    pval_mask = de_results["pvals_adj"] < pval_thresh

    if only_positive:
        fc_mask = de_results["logfoldchanges"] > min_log_fc
    else:
        fc_mask = de_results["logfoldchanges"].abs() > min_log_fc

    filtered = de_results.loc[pval_mask & fc_mask].copy()
    filtered.sort_values("logfoldchanges", ascending=False, inplace=True)

    filtered["z"] = filtered["logfoldchanges"] * (
        -np.log10(filtered["pvals_adj"] + 1e-9)
    )

    return filtered


def _run_linear_model_de(
    adata_subset: AnnData,
    groupby_variable: str,
    group_selection: str,
    covariates: Sequence[str],
    use_raw: bool = False,
    gene_chunk_size: int = _OLS_GENE_CHUNK,
) -> pd.DataFrame:
    obs = adata_subset.obs

    is_group = (obs[groupby_variable] == group_selection).astype(float).values

    dummy_blocks: list[np.ndarray] = []
    cov_sizes: list[int] = []
    for cov in covariates:
        dummies = pd.get_dummies(
            obs[cov].astype(str), prefix=cov, drop_first=True, dtype=float
        )
        dummy_blocks.append(dummies.values)
        cov_sizes.append(dummies.shape[1])

    pieces = [np.ones((len(obs), 1)), is_group[:, None]]
    for block in dummy_blocks:
        pieces.append(block)
    X = np.column_stack(pieces)
    del dummy_blocks, pieces

    covariate_desc = ", ".join(
        f"{sz} {c}" for c, sz in zip(covariates, cov_sizes)
    )
    logger.info(
        "Design matrix: %d cells × %d covariates "
        "(1 intercept + 1 group + %s)",
        X.shape[0], X.shape[1], covariate_desc,
    )

    XtX_inv = np.linalg.pinv(X.T @ X)
    n, k = X.shape

    source = adata_subset.raw if use_raw else adata_subset
    expr_mat = source.X
    is_sparse = sp_sparse.issparse(expr_mat)
    n_genes = expr_mat.shape[1]
    gene_names = source.var_names

    all_betas = np.empty(n_genes, dtype=np.float64)
    all_t = np.empty(n_genes, dtype=np.float64)
    all_pvals = np.empty(n_genes, dtype=np.float64)

    for start in range(0, n_genes, gene_chunk_size):
        end = min(start + gene_chunk_size, n_genes)
        Y_chunk = expr_mat[:, start:end]
        if is_sparse:
            Y_chunk = Y_chunk.toarray()
        elif not isinstance(Y_chunk, np.ndarray):
            Y_chunk = np.asarray(Y_chunk)

        beta_chunk = XtX_inv @ (X.T @ Y_chunk)
        residuals = Y_chunk - X @ beta_chunk
        sigma2 = (residuals ** 2).sum(axis=0) / (n - k)

        se_group = np.sqrt(XtX_inv[1, 1] * sigma2)
        with np.errstate(divide="ignore", invalid="ignore"):
            t_stat = beta_chunk[1, :] / se_group
        t_stat = np.nan_to_num(t_stat, nan=0.0)
        pvals = 2 * scipy_stats.t.sf(np.abs(t_stat), df=n - k)
        pvals = np.nan_to_num(pvals, nan=1.0)

        all_betas[start:end] = beta_chunk[1, :]
        all_t[start:end] = t_stat
        all_pvals[start:end] = pvals

        del Y_chunk, beta_chunk, residuals, sigma2, se_group, t_stat, pvals

    _, pvals_adj, _, _ = multipletests(all_pvals, method="fdr_bh")

    de_results = pd.DataFrame({
        "names": gene_names,
        "coef": all_betas,
        "logfoldchanges": all_betas / np.log(2),
        "t_statistic": all_t,
        "pvals": all_pvals,
        "pvals_adj": pvals_adj,
    })

    return de_results


def _lodo_single_fold_from_indices(
    adata: AnnData,
    keep_indices: np.ndarray,
    level: str,
    groupby_variable: str,
    group_selection: str,
    reference: Optional[str],
    method: str,
    pval_thresh: float,
    min_log_fc: float,
    only_positive: bool,
    use_raw: bool,
    min_pct_nz_group: float = 0.0,
    min_pct_nz_difference: float = 0.0,
    pts: bool = True,
) -> Optional[pd.DataFrame]:
    adata_fold = adata[keep_indices].copy()

    de_kwargs = dict(
        groupby=groupby_variable,
        groups=[group_selection],
        method=method,
        use_raw=use_raw,
        pts=pts,
    )
    if reference is not None:
        de_kwargs["reference"] = reference

    sc.tl.rank_genes_groups(adata_fold, **de_kwargs)
    iter_results = sc.get.rank_genes_groups_df(adata_fold, group=group_selection)

    del adata_fold
    gc.collect()

    pval_mask = iter_results["pvals_adj"] < pval_thresh
    if only_positive:
        lfc_mask = iter_results["logfoldchanges"] > min_log_fc
    else:
        lfc_mask = iter_results["logfoldchanges"].abs() > min_log_fc
    combined_mask = pval_mask & lfc_mask
    if "pct_nz_group" in iter_results.columns:
        combined_mask = combined_mask & (iter_results["pct_nz_group"] > min_pct_nz_group)
    iter_sig = iter_results[combined_mask].copy()

    if min_pct_nz_difference > 0.0 and "pct_nz_reference" in iter_sig.columns:
        pct_nz_diff = iter_sig["pct_nz_group"] - iter_sig["pct_nz_reference"]
        drop_mask = (iter_sig["pct_nz_reference"] < 0.01) & (
            pct_nz_diff < min_pct_nz_difference
        )
        iter_sig = iter_sig.loc[~drop_mask]

    iter_sig["_excluded_level"] = level
    del iter_results
    return iter_sig


def leave_one_out_de(
    adata: AnnData,
    leave_out_variable: str,
    groupby_variable: str,
    group_selection: str,
    *,
    reference: Optional[str] = None,
    method: str = "wilcoxon",
    min_cell_count: int = 5,
    pval_thresh: float = 0.01,
    min_log_fc: float = 0.25,
    only_positive: bool = True,
    robustness_thresh: float = 0.5,
    use_raw: bool = False,
    min_pct_nz_group: float = 0.0,
    min_pct_nz_difference: float = 0.0,
    pts: bool = True,
    n_jobs: int = 1,
) -> pd.DataFrame:
    levels = sorted(adata.obs[leave_out_variable].unique())
    n_workers = _resolve_n_jobs(n_jobs)

    logger.info(
        "Leave-one-out DE: %d levels of '%s' (n_jobs=%d → %d workers)",
        len(levels), leave_out_variable, n_jobs, n_workers,
    )

    obs_var = adata.obs[leave_out_variable]
    obs_grp = adata.obs[groupby_variable]
    fold_specs: list[tuple[str, np.ndarray]] = []

    for level in levels:
        keep_mask = obs_var != level
        keep_idx = np.where(keep_mask.values)[0]

        grp_vals = obs_grp.iloc[keep_idx]
        group_counts = grp_vals.value_counts()
        n_test = group_counts.get(group_selection, 0)
        n_ref = group_counts.drop(group_selection, errors="ignore").sum()

        if n_test < min_cell_count or n_ref < min_cell_count:
            logger.info(
                "  Excluding %s: SKIPPED (test=%d, ref=%d)",
                level, n_test, n_ref,
            )
            continue

        logger.info("  Excluding %s: test=%d, ref=%d", level, n_test, n_ref)
        fold_specs.append((level, keep_idx))

    n_valid = len(fold_specs)
    logger.info("Valid folds: %d / %d", n_valid, len(levels))

    if n_valid == 0:
        logger.warning("No valid iterations — returning empty DataFrame")
        return pd.DataFrame(
            columns=[
                "names", "n_sig", "n_iterations", "frac_sig",
                "direction_consistent", "direction_up",
            ]
        )

    def _run_fold(level: str, keep_idx: np.ndarray) -> Optional[pd.DataFrame]:
        return _lodo_single_fold_from_indices(
            adata, keep_idx, level,
            groupby_variable, group_selection, reference,
            method, pval_thresh, min_log_fc, only_positive, use_raw,
            min_pct_nz_group, min_pct_nz_difference, pts,
        )

    if n_workers > 1:
        fold_results = Parallel(n_jobs=n_workers, backend="loky")(
            delayed(_run_fold)(level, keep_idx)
            for level, keep_idx in fold_specs
        )
    else:
        fold_results = [_run_fold(level, keep_idx) for level, keep_idx in fold_specs]

    valid_results = [r for r in fold_results if r is not None]
    del fold_results

    if not valid_results:
        logger.warning("No valid iterations — returning empty DataFrame")
        return pd.DataFrame(
            columns=[
                "names", "n_sig", "n_iterations", "frac_sig",
                "direction_consistent", "direction_up",
            ]
        )

    sig_counts: Dict[str, int] = defaultdict(int)
    directions: Dict[str, List[bool]] = defaultdict(list)

    for iter_sig in valid_results:
        names = iter_sig["names"].values
        lfc_values = iter_sig["logfoldchanges"].values
        is_up = lfc_values > 0
        for gene, up in zip(names, is_up):
            sig_counts[gene] += 1
            directions[gene].append(bool(up))

    del valid_results
    gc.collect()

    if n_valid <= 1:
        logger.warning(
            "Only %d valid LODO fold(s) — cannot compute robustness; "
            "returning empty DataFrame",
            n_valid,
        )
        return pd.DataFrame(
            columns=[
                "names", "n_sig", "n_iterations", "frac_sig",
                "direction_consistent", "direction_up",
            ]
        )

    rows = []
    for gene in sig_counts:
        n_sig = sig_counts[gene]
        dirs = directions[gene]
        consistent = all(dirs) or not any(dirs)
        up = all(dirs) if consistent else None
        rows.append({
            "names": gene,
            "n_sig": n_sig,
            "n_iterations": n_valid,
            "frac_sig": n_sig / n_valid,
            "direction_consistent": consistent,
            "direction_up": up,
        })

    lodo_df = pd.DataFrame(
        rows,
        columns=[
            "names", "n_sig", "n_iterations", "frac_sig",
            "direction_consistent", "direction_up",
        ],
    )
    if len(lodo_df) > 0:
        lodo_df = lodo_df.sort_values("frac_sig", ascending=False)
    return lodo_df


def _per_donor_single_fold_from_indices(
    adata: AnnData,
    donor_indices: np.ndarray,
    donor_label: str,
    groupby_variable: str,
    group_selection: str,
    reference: Optional[str],
    method: str,
    pval_thresh: float,
    min_log_fc: float,
    only_positive: bool,
    use_raw: bool,
    min_pct_nz_group: float = 0.0,
    min_pct_nz_difference: float = 0.0,
    pts: bool = True,
) -> Optional[pd.DataFrame]:
    adata_fold = adata[donor_indices].copy()

    de_kwargs = dict(
        groupby=groupby_variable,
        groups=[group_selection],
        method=method,
        use_raw=use_raw,
        pts=pts,
    )
    if reference is not None:
        de_kwargs["reference"] = reference

    sc.tl.rank_genes_groups(adata_fold, **de_kwargs)
    iter_results = sc.get.rank_genes_groups_df(adata_fold, group=group_selection)

    del adata_fold
    gc.collect()

    pval_mask = iter_results["pvals_adj"] < pval_thresh
    if only_positive:
        lfc_mask = iter_results["logfoldchanges"] > min_log_fc
    else:
        lfc_mask = iter_results["logfoldchanges"].abs() > min_log_fc
    combined_mask = pval_mask & lfc_mask
    if "pct_nz_group" in iter_results.columns:
        combined_mask = combined_mask & (iter_results["pct_nz_group"] > min_pct_nz_group)
    iter_sig = iter_results[combined_mask].copy()

    if min_pct_nz_difference > 0.0 and "pct_nz_reference" in iter_sig.columns:
        pct_nz_diff = iter_sig["pct_nz_group"] - iter_sig["pct_nz_reference"]
        drop_mask = (iter_sig["pct_nz_reference"] < 0.01) & (
            pct_nz_diff < min_pct_nz_difference
        )
        iter_sig = iter_sig.loc[~drop_mask]

    iter_sig["_donor"] = donor_label
    del iter_results
    return iter_sig


def per_donor_de(
    adata: AnnData,
    donor_variable: str,
    groupby_variable: str,
    group_selection: str,
    *,
    reference: Optional[str] = None,
    method: str = "wilcoxon",
    min_cell_count: int = 5,
    pval_thresh: float = 0.01,
    min_log_fc: float = 0.25,
    only_positive: bool = True,
    consistency_thresh: float = 0.5,
    use_raw: bool = False,
    min_pct_nz_group: float = 0.0,
    min_pct_nz_difference: float = 0.0,
    pts: bool = True,
    n_jobs: int = 1,
) -> pd.DataFrame:
    donors = sorted(adata.obs[donor_variable].unique())
    n_workers = _resolve_n_jobs(n_jobs)

    logger.info(
        "Per-donor DE: %d donors in '%s' (n_jobs=%d → %d workers)",
        len(donors), donor_variable, n_jobs, n_workers,
    )

    obs_donor = adata.obs[donor_variable]
    obs_grp = adata.obs[groupby_variable]
    fold_specs: list[tuple[str, np.ndarray]] = []

    for donor in donors:
        donor_mask = obs_donor == donor
        donor_idx = np.where(donor_mask.values)[0]

        grp_vals = obs_grp.iloc[donor_idx]
        group_counts = grp_vals.value_counts()
        n_test = group_counts.get(group_selection, 0)
        n_ref = group_counts.drop(group_selection, errors="ignore").sum()

        if n_test < min_cell_count or n_ref < min_cell_count:
            logger.info(
                "  Donor %s: SKIPPED (test=%d, ref=%d)",
                donor, n_test, n_ref,
            )
            continue

        logger.info("  Donor %s: test=%d, ref=%d", donor, n_test, n_ref)
        fold_specs.append((donor, donor_idx))

    n_valid = len(fold_specs)
    logger.info("Valid donors: %d / %d", n_valid, len(donors))

    if n_valid == 0:
        logger.warning("No valid donors — returning empty DataFrame")
        return pd.DataFrame(
            columns=[
                "names", "n_sig", "n_donors", "frac_sig",
                "direction_consistent", "direction_up",
            ]
        )

    def _run_fold(donor: str, donor_idx: np.ndarray) -> Optional[pd.DataFrame]:
        return _per_donor_single_fold_from_indices(
            adata, donor_idx, donor,
            groupby_variable, group_selection, reference,
            method, pval_thresh, min_log_fc, only_positive, use_raw,
            min_pct_nz_group, min_pct_nz_difference, pts,
        )

    if n_workers > 1:
        fold_results = Parallel(n_jobs=n_workers, backend="loky")(
            delayed(_run_fold)(donor, donor_idx)
            for donor, donor_idx in fold_specs
        )
    else:
        fold_results = [_run_fold(donor, donor_idx) for donor, donor_idx in fold_specs]

    valid_results = [r for r in fold_results if r is not None]
    del fold_results

    if not valid_results:
        logger.warning("No valid donors — returning empty DataFrame")
        return pd.DataFrame(
            columns=[
                "names", "n_sig", "n_donors", "frac_sig",
                "direction_consistent", "direction_up",
            ]
        )

    sig_counts: Dict[str, int] = defaultdict(int)
    directions: Dict[str, List[bool]] = defaultdict(list)

    for iter_sig in valid_results:
        names = iter_sig["names"].values
        lfc_values = iter_sig["logfoldchanges"].values
        is_up = lfc_values > 0
        for gene, up in zip(names, is_up):
            sig_counts[gene] += 1
            directions[gene].append(bool(up))

    del valid_results
    gc.collect()

    if n_valid <= 1:
        logger.warning(
            "Only %d valid donor(s) — cannot compute consistency; "
            "returning empty DataFrame",
            n_valid,
        )
        return pd.DataFrame(
            columns=[
                "names", "n_sig", "n_donors", "frac_sig",
                "direction_consistent", "direction_up",
            ]
        )

    rows = []
    for gene in sig_counts:
        n_sig = sig_counts[gene]
        dirs = directions[gene]
        consistent = all(dirs) or not any(dirs)
        up = all(dirs) if consistent else None
        rows.append({
            "names": gene,
            "n_sig": n_sig,
            "n_donors": n_valid,  
            "frac_sig": n_sig / n_valid,
            "direction_consistent": consistent,
            "direction_up": up,
        })

    per_donor_df = pd.DataFrame(
        rows,
        columns=[
            "names", "n_sig", "n_donors", "frac_sig",
            "direction_consistent", "direction_up",
        ],
    )
    if len(per_donor_df) > 0:
        per_donor_df = per_donor_df.sort_values("frac_sig", ascending=False)
    return per_donor_df


def find_senescence_associated_genes(
    adata: AnnData,
    variable_key: str,
    variable_value: str,
    groupby_variable: str,
    *,
    group_selection: str = "Positive",
    covariates: Optional[Sequence[str]] = None,
    leave_one_out_variable: Optional[str] = None,
    lodo_reference: Optional[str] = None,
    lodo_robustness_thresh: float = 0.5,
    lodo_pval_thresh: Optional[float] = None,
    lodo_n_jobs: int = 1,
    per_donor_variable: Optional[str] = None,
    per_donor_consistency_thresh: float = 0.5,
    per_donor_pval_thresh: Optional[float] = None,
    per_donor_n_jobs: int = 1,
    key_added: str = "rank_gene_groups",
    method: str = "wilcoxon",
    min_cell_count: int = 5,
    only_positive: bool = True,
    min_log_fc: float = 0.25,
    pval_thresh: float = 0.01,
    min_pct_nz_group: float = 0.0,
    min_pct_nz_difference: float = 0.0,
    return_format: Literal["names", "z", "lfc", "df"] = "names",
    pts: bool = True,
    use_raw: bool = False,
    output_path: Optional[Union[str, Path]] = None,
    subsample_negatives_to: Optional[int] = None,
    subsample_targets: Optional[Dict[str, int]] = None,
    random_state: int = 42,
) -> MarkerResult:
    subset_mask = adata.obs[variable_key] == variable_value
    adata_subset = adata[subset_mask]
    logger.info(
        "%s: %s  —  shape %s", variable_key, variable_value, adata_subset.shape
    )

    # Determine subsampling target for this category
    _subsample_n = subsample_negatives_to
    if subsample_targets is not None and variable_value in subsample_targets:
        _subsample_n = subsample_targets[variable_value]

    # Downsample non-senescent ("Negative") cells if requested
    if _subsample_n is not None:
        neg_mask = adata_subset.obs[groupby_variable] != group_selection
        pos_mask = adata_subset.obs[groupby_variable] == group_selection
        n_neg = int(neg_mask.sum())
        if n_neg > _subsample_n:
            neg_idx = np.where(neg_mask.values)[0]
            rng = np.random.RandomState(random_state)
            keep_neg = rng.choice(neg_idx, size=_subsample_n, replace=False)
            pos_idx = np.where(pos_mask.values)[0]
            keep_all = np.sort(np.concatenate([pos_idx, keep_neg]))
            adata_subset = adata_subset[keep_all].copy()
            logger.info(
                "  Subsampled negatives: %d → %d (positives: %d)",
                n_neg, _subsample_n, len(pos_idx),
            )

    group_counts = adata_subset.obs[groupby_variable].value_counts()
    if len(group_counts) < 2 or group_counts.min() < min_cell_count:
        logger.warning(
            "Skipping %s=%s (group counts: %s)",
            variable_key, variable_value, group_counts.to_dict(),
        )
        return None

    if covariates is not None and len(covariates) > 0:
        lm_results_full = _run_linear_model_de(
            adata_subset,
            groupby_variable=groupby_variable,
            group_selection=group_selection,
            covariates=covariates,
            use_raw=use_raw,
        )
        markers = filter_lm_results(
            lm_results_full,
            only_positive=only_positive,
            pval_thresh=pval_thresh,
            min_log_fc=min_log_fc,
        )
        del lm_results_full
        gc.collect()

        if leave_one_out_variable is not None:
            lodo_pval = lodo_pval_thresh if lodo_pval_thresh is not None else pval_thresh

            lodo_df = leave_one_out_de(
                adata_subset,
                leave_out_variable=leave_one_out_variable,
                groupby_variable=groupby_variable,
                group_selection=group_selection,
                reference=lodo_reference,
                method=method,
                min_cell_count=min_cell_count,
                pval_thresh=lodo_pval,
                min_log_fc=min_log_fc,
                only_positive=only_positive,
                robustness_thresh=lodo_robustness_thresh,
                use_raw=use_raw,
                min_pct_nz_group=min_pct_nz_group,
                min_pct_nz_difference=min_pct_nz_difference,
                pts=pts,
                n_jobs=lodo_n_jobs,
            )

            adata_wx = adata_subset.copy()
            wilcoxon_key = f"{key_added}_{groupby_variable}_wilcoxon"
            sc.tl.rank_genes_groups(
                adata_wx,
                groupby=groupby_variable,
                key_added=wilcoxon_key,
                use_raw=use_raw,
                pts=pts,
                method=method,
            )
            wilcoxon_markers = sc.get.rank_genes_groups_df(
                adata_wx,
                key=wilcoxon_key,
                group=group_selection,
                pval_cutoff=pval_thresh,
            )
            del adata_wx
            gc.collect()

            wilcoxon_markers = filter_markers(
                wilcoxon_markers,
                only_positive=only_positive,
                pval_thresh=pval_thresh,
                min_log_fc=min_log_fc,
                min_pct_nz_group=min_pct_nz_group,
                min_pct_nz_difference=min_pct_nz_difference,
            )

            per_donor_df = None
            if per_donor_variable is not None:
                pd_pval = (
                    per_donor_pval_thresh
                    if per_donor_pval_thresh is not None
                    else pval_thresh
                )
                per_donor_df = per_donor_de(
                    adata_subset,
                    donor_variable=per_donor_variable,
                    groupby_variable=groupby_variable,
                    group_selection=group_selection,
                    reference=lodo_reference,
                    method=method,
                    min_cell_count=min_cell_count,
                    pval_thresh=pd_pval,
                    min_log_fc=min_log_fc,
                    only_positive=only_positive,
                    consistency_thresh=per_donor_consistency_thresh,
                    use_raw=use_raw,
                    min_pct_nz_group=min_pct_nz_group,
                    min_pct_nz_difference=min_pct_nz_difference,
                    pts=pts,
                    n_jobs=per_donor_n_jobs,
                )

            logger.info(
                "LM significant: %d | Wilcoxon (no covariates): %d",
                len(markers), len(wilcoxon_markers),
            )

            result_dict = {
                "lm": markers,
                "lodo": lodo_df,
                "wilcoxon": wilcoxon_markers,
            }
            if per_donor_df is not None:
                result_dict["per_donor"] = per_donor_df

            if output_path is not None:
                out_dir = Path(output_path)
                out_dir.mkdir(parents=True, exist_ok=True)
                markers.to_csv(
                    out_dir / f"{variable_value}_lm_markers.csv", index=False
                )
                if len(lodo_df) > 0:
                    lodo_df.to_csv(
                        out_dir / f"{variable_value}_lodo_summary.csv", index=False
                    )
                else:
                    logger.warning(
                        "%s: no valid LODO iterations — skipping LODO CSV",
                        variable_value,
                    )
                wilcoxon_markers.to_csv(
                    out_dir / f"{variable_value}_wilcoxon_markers.csv", index=False
                )
                if per_donor_df is not None:
                    if len(per_donor_df) > 0:
                        per_donor_df.to_csv(
                            out_dir / f"{variable_value}_per_donor_summary.csv",
                            index=False,
                        )
                    else:
                        logger.warning(
                            "%s: no valid per-donor DE iterations — skipping "
                            "per-donor CSV",
                            variable_value,
                        )

                has_data = any(
                    isinstance(v, pd.DataFrame) and len(v) > 0
                    for v in result_dict.values()
                )
                return result_dict if has_data else None

            return result_dict

        if output_path is not None:
            out_dir = Path(output_path)
            out_dir.mkdir(parents=True, exist_ok=True)
            csv_path = out_dir / f"{variable_value}_lm_markers.csv"
            markers.to_csv(csv_path, index=False)
            logger.info("Saved LM markers → %s", csv_path)
            return markers if len(markers) > 0 else None

        return _format_results(markers, return_format)

    else:
        adata_de = adata_subset.copy()
        del adata_subset

        result_key = f"{key_added}_{groupby_variable}"
        sc.tl.rank_genes_groups(
            adata_de,
            groupby=groupby_variable,
            key_added=result_key,
            use_raw=use_raw,
            pts=pts,
            method=method,
        )
        markers = sc.get.rank_genes_groups_df(
            adata_de,
            key=result_key,
            group=group_selection,
            pval_cutoff=pval_thresh,
        )

        del adata_de
        gc.collect()

        markers = filter_markers(
            markers,
            only_positive=only_positive,
            pval_thresh=pval_thresh,
            min_log_fc=min_log_fc,
            min_pct_nz_group=min_pct_nz_group,
            min_pct_nz_difference=min_pct_nz_difference,
        )

    if output_path is not None:
        out_dir = Path(output_path)
        out_dir.mkdir(parents=True, exist_ok=True)
        csv_path = out_dir / f"{variable_value}_senescence_markers.csv"
        markers.to_csv(csv_path, index=False)
        logger.info("Saved markers → %s", csv_path)
        return markers if len(markers) > 0 else None

    return _format_results(markers, return_format)


def _diagnose_skip_reason(
    adata: AnnData,
    variable_key: str,
    variable_value: str,
    groupby_variable: str,
    group_selection: str,
    min_cell_count: int,
    result: MarkerResult,
    output_path: Optional[Union[str, Path]],
) -> Optional[str]:
    if result is None:
        subset_mask = adata.obs[variable_key] == variable_value
        n_cells = int(subset_mask.sum())
        if n_cells == 0:
            return f"no cells found for {variable_key}={variable_value}"
        group_counts = adata.obs.loc[subset_mask, groupby_variable].value_counts()
        if len(group_counts) < 2:
            return (
                f"fewer than 2 groups in '{groupby_variable}' "
                f"(found: {group_counts.to_dict()})"
            )
        if group_selection not in group_counts.index:
            return (
                f"group '{group_selection}' not present "
                f"(found: {group_counts.to_dict()})"
            )
        if group_counts.min() < min_cell_count:
            return (
                f"a group has < {min_cell_count} cells "
                f"(counts: {group_counts.to_dict()})"
            )
        return "no genes passed filtering thresholds"

    if isinstance(result, pd.DataFrame) and len(result) == 0:
        return "no genes passed filtering thresholds"
    if isinstance(result, list) and len(result) == 0:
        return "no genes passed filtering thresholds"
    if isinstance(result, dict):
        if not any(
            (isinstance(v, pd.DataFrame) and len(v) > 0)
            or (isinstance(v, list) and len(v) > 0)
            for v in result.values()
        ):
            parts = [
                f"{k}: 0 genes"
                for k, v in result.items()
                if isinstance(v, (pd.DataFrame, list))
            ]
            return (
                "no genes passed filtering in any analysis mode "
                f"({', '.join(parts)})"
            )

    return None


def _process_single_category(
    adata: AnnData,
    variable: str,
    variable_key: str,
    kwargs: dict,
) -> Tuple[str, MarkerResult, Optional[str]]:
    result = find_senescence_associated_genes(
        adata,
        variable_key=variable_key,
        variable_value=variable,
        **kwargs,
    )
    skip_reason = _diagnose_skip_reason(
        adata, variable_key, variable,
        kwargs["groupby_variable"],
        kwargs["group_selection"],
        kwargs["min_cell_count"],
        result,
        kwargs.get("output_path"),
    )
    return variable, result, skip_reason


def identify_robust_sags(
    results: Dict[str, dict],
    *,
    lodo_frac_thresh: float = 0.5,
    per_donor_frac_thresh: float = 0.5,
    min_methods: int = 2,
    exclude_methods: Optional[List[str]] = None,
    output_path: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    rows = []

    for tissue_ct, result_dict in results.items():
        if not isinstance(result_dict, dict):
            continue

        wx_df = result_dict.get("wilcoxon")
        if wx_df is None or not isinstance(wx_df, pd.DataFrame) or len(wx_df) == 0:
            continue

        wx_genes = set(wx_df["names"].values)
        wx_lookup = wx_df.set_index("names")

        lm_df = result_dict.get("lm")
        lm_genes = set()
        if lm_df is not None and isinstance(lm_df, pd.DataFrame) and len(lm_df) > 0:
            lm_genes = set(lm_df["names"].values)

        lodo_df = result_dict.get("lodo")
        lodo_genes = set()
        lodo_lookup: Union[pd.DataFrame, dict] = {}
        if lodo_df is not None and isinstance(lodo_df, pd.DataFrame) and len(lodo_df) > 0:
            lodo_pass = lodo_df[
                (lodo_df["n_iterations"] > 1)
                & (lodo_df["frac_sig"] >= lodo_frac_thresh)
                & (lodo_df["direction_consistent"] == True)
                & (lodo_df["direction_up"] == True)
            ]
            lodo_genes = set(lodo_pass["names"].values)
            lodo_lookup = lodo_df.set_index("names")

        pd_df = result_dict.get("per_donor")
        pd_genes = set()
        pd_lookup: Union[pd.DataFrame, dict] = {}
        if pd_df is not None and isinstance(pd_df, pd.DataFrame) and len(pd_df) > 0:
            pd_pass = pd_df[
                (pd_df["n_donors"] > 1)
                & (pd_df["frac_sig"] >= per_donor_frac_thresh)
                & (pd_df["direction_consistent"] == True)
                & (pd_df["direction_up"] == True)
            ]
            pd_genes = set(pd_pass["names"].values)
            pd_lookup = pd_df.set_index("names")

        _excluded = set(exclude_methods or [])

        for gene in wx_genes:
            in_lm = gene in lm_genes
            in_lodo = gene in lodo_genes
            in_pd = gene in pd_genes
            n_methods = 1  # wilcoxon always counts
            if "lm" not in _excluded:
                n_methods += int(in_lm)
            if "lodo" not in _excluded:
                n_methods += int(in_lodo)
            if "per_donor" not in _excluded:
                n_methods += int(in_pd)

            if n_methods < min_methods:
                continue
                
            lodo_frac = np.nan
            lodo_n_iter = np.nan
            if isinstance(lodo_lookup, pd.DataFrame) and gene in lodo_lookup.index:
                lodo_frac = lodo_lookup.loc[gene, "frac_sig"]
                lodo_n_iter = lodo_lookup.loc[gene, "n_iterations"]

            pd_frac = np.nan
            pd_n_donors = np.nan
            if isinstance(pd_lookup, pd.DataFrame) and gene in pd_lookup.index:
                pd_frac = pd_lookup.loc[gene, "frac_sig"]
                pd_n_donors = pd_lookup.loc[gene, "n_donors"]

            wx_logfc = wx_lookup.loc[gene, "logfoldchanges"]
            wx_pval = wx_lookup.loc[gene, "pvals_adj"]

            rows.append({
                "tissue_celltype": tissue_ct,
                "gene": gene,
                "in_wilcoxon": True,
                "in_lm": in_lm,
                "in_lodo": in_lodo,
                "in_per_donor": in_pd,
                "n_methods": n_methods,
                "lodo_frac_sig": lodo_frac,
                "lodo_n_iterations": lodo_n_iter,
                "per_donor_frac_sig": pd_frac,
                "per_donor_n_donors": pd_n_donors,
                "wilcoxon_logfc": wx_logfc,
                "wilcoxon_pval_adj": wx_pval,
            })

    robust_df = pd.DataFrame(
        rows,
        columns=[
            "tissue_celltype", "gene", "in_wilcoxon", "in_lm", "in_lodo",
            "in_per_donor", "n_methods", "lodo_frac_sig", "lodo_n_iterations",
            "per_donor_frac_sig", "per_donor_n_donors", "wilcoxon_logfc",
            "wilcoxon_pval_adj",
        ],
    )

    if len(robust_df) > 0:
        robust_df.sort_values(
            ["n_methods", "wilcoxon_logfc"],
            ascending=[False, False],
            inplace=True,
        )

    n_genes = robust_df["gene"].nunique() if len(robust_df) > 0 else 0
    n_categories = robust_df["tissue_celltype"].nunique() if len(robust_df) > 0 else 0
    logger.info(
        "Robust SAGs: %d unique genes across %d tissue-celltypes "
        "(%d gene–tissue pairs)",
        n_genes, n_categories, len(robust_df),
    )

    if output_path is not None:
        out_dir = Path(output_path)
        out_dir.mkdir(parents=True, exist_ok=True)
        csv_path = out_dir / "robust_sags.csv"
        robust_df.to_csv(csv_path, index=False)
        logger.info("Saved robust SAGs → %s", csv_path)

    return robust_df


def summarize_robust_sags(
    robust_df: pd.DataFrame,
    *,
    output_path: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    if len(robust_df) == 0:
        logger.warning("No robust SAGs to summarize")
        return pd.DataFrame(
            columns=[
                "gene", "n_tissue_celltypes", "tissue_celltypes",
                "mean_n_methods", "mean_wilcoxon_logfc",
            ]
        )

    summary = (
        robust_df.groupby("gene")
        .agg(
            n_tissue_celltypes=("tissue_celltype", "nunique"),
            tissue_celltypes=("tissue_celltype", lambda x: ", ".join(sorted(x.unique()))),
            mean_n_methods=("n_methods", "mean"),
            mean_wilcoxon_logfc=("wilcoxon_logfc", "mean"),
        )
        .reset_index()
        .sort_values("n_tissue_celltypes", ascending=False)
    )

    logger.info(
        "SAG summary: %d unique genes, top gene appears in %d tissue-celltypes",
        len(summary),
        summary["n_tissue_celltypes"].iloc[0] if len(summary) > 0 else 0,
    )

    if output_path is not None:
        out_dir = Path(output_path)
        out_dir.mkdir(parents=True, exist_ok=True)
        csv_path = out_dir / "robust_sags_summary.csv"
        summary.to_csv(csv_path, index=False)
        logger.info("Saved SAG summary → %s", csv_path)

    return summary


def process_adata_for_sags(
    adata: AnnData,
    variable_key: str,
    *,
    groupby_variable: str = "CDKN2A+ MKI67-",
    group_selection: str = "Positive",
    covariates: Optional[Sequence[str]] = None,
    leave_one_out_variable: Optional[str] = None,
    lodo_reference: Optional[str] = None,
    lodo_robustness_thresh: float = 0.5,
    lodo_pval_thresh: Optional[float] = None,
    lodo_n_jobs: int = 1,
    per_donor_variable: Optional[str] = None,
    per_donor_consistency_thresh: float = 0.5,
    per_donor_pval_thresh: Optional[float] = None,
    per_donor_n_jobs: int = 1,
    key_added: str = "rank_gene_groups",
    method: str = "wilcoxon",
    pts: bool = True,
    only_positive: bool = True,
    min_log_fc: float = 0.25,
    pval_thresh: float = 0.01,
    min_cell_count: int = 5,
    min_pct_nz_group: float = 0.0,
    min_pct_nz_difference: float = 0.0,
    use_raw: bool = False,
    return_format: Literal["names", "z", "lfc", "df"] = "names",
    output_path: Optional[Union[str, Path]] = None,
    n_jobs: int = 1,
    identify_robust: bool = False,
    robust_lodo_frac_thresh: float = 0.5,
    robust_per_donor_frac_thresh: float = 0.5,
    robust_min_methods: int = 2,
    robust_exclude_methods: Optional[List[str]] = None,
    subsample_negatives_to: Optional[int] = None,
    subsample_targets: Optional[Dict[str, int]] = None,
    random_state: int = 42,
) -> Dict[str, MarkerResult]:
    col = adata.obs[variable_key]
    categories = list(col.cat.categories if hasattr(col, "cat") else sorted(col.unique()))
    n_workers = _resolve_n_jobs(n_jobs)

    logger.info(
        "Processing %d categories in '%s' (n_jobs=%d → %d workers)",
        len(categories), variable_key, n_jobs, n_workers,
    )

    shared_kwargs = dict(
        groupby_variable=groupby_variable,
        group_selection=group_selection,
        covariates=covariates,
        leave_one_out_variable=leave_one_out_variable,
        lodo_reference=lodo_reference,
        lodo_robustness_thresh=lodo_robustness_thresh,
        lodo_pval_thresh=lodo_pval_thresh,
        lodo_n_jobs=lodo_n_jobs,
        per_donor_variable=per_donor_variable,
        per_donor_consistency_thresh=per_donor_consistency_thresh,
        per_donor_pval_thresh=per_donor_pval_thresh,
        per_donor_n_jobs=per_donor_n_jobs,
        key_added=key_added,
        method=method,
        min_cell_count=min_cell_count,
        only_positive=only_positive,
        min_log_fc=min_log_fc,
        pval_thresh=pval_thresh,
        min_pct_nz_group=min_pct_nz_group,
        min_pct_nz_difference=min_pct_nz_difference,
        return_format=return_format,
        pts=pts,
        use_raw=use_raw,
        output_path=output_path,
        subsample_negatives_to=subsample_negatives_to,
        subsample_targets=subsample_targets,
        random_state=random_state,
    )

    if n_workers > 1:
        results = Parallel(n_jobs=n_workers, backend="loky")(
            delayed(_process_single_category)(
                adata, variable, variable_key, shared_kwargs,
            )
            for variable in categories
        )
    else:
        results = [
            _process_single_category(adata, variable, variable_key, shared_kwargs)
            for variable in categories
        ]

    marker_dict: Dict[str, MarkerResult] = {}
    saved_categories: List[str] = []
    skipped_categories: List[Tuple[str, str]] = []

    for variable, result, skip_reason in results:
        if skip_reason is not None:
            skipped_categories.append((variable, skip_reason))
        if result is None:
            if output_path is not None and skip_reason is None:
                saved_categories.append(variable)
            continue
        if output_path is not None:
            saved_categories.append(variable)
            continue
        if isinstance(result, dict):
            if any(
                (isinstance(v, pd.DataFrame) and len(v) > 0)
                or (isinstance(v, list) and len(v) > 0)
                for v in result.values()
            ):
                marker_dict[variable] = result
        elif isinstance(result, pd.DataFrame) and len(result) > 0:
            marker_dict[variable] = result
        elif isinstance(result, list) and len(result) > 0:
            marker_dict[variable] = result

    del results
    gc.collect()

    if skipped_categories:
        logger.warning(
            "No DGEA results for %d / %d categories:",
            len(skipped_categories), len(categories),
        )
        for cat, reason in skipped_categories:
            logger.warning("  %s=%s — %s", variable_key, cat, reason)

    if output_path is not None:
        logger.info(
            "Saved results for %d / %d categories to %s",
            len(saved_categories), len(categories), output_path,
        )
        if not identify_robust:
            return {}
        out_dir = Path(output_path)
        for cat in saved_categories:
            cat_dict = {}
            for suffix, key in [
                ("_wilcoxon_markers.csv", "wilcoxon"),
                ("_lm_markers.csv", "lm"),
                ("_lodo_summary.csv", "lodo"),
                ("_per_donor_summary.csv", "per_donor"),
            ]:
                csv_path = out_dir / f"{cat}{suffix}"
                if csv_path.exists():
                    cat_dict[key] = pd.read_csv(csv_path)
            if cat_dict:
                marker_dict[cat] = cat_dict

    logger.info(
        "Found markers for %d / %d categories", len(marker_dict), len(categories)
    )

    if identify_robust and marker_dict:
        robust_df = identify_robust_sags(
            marker_dict,
            lodo_frac_thresh=robust_lodo_frac_thresh,
            per_donor_frac_thresh=robust_per_donor_frac_thresh,
            min_methods=robust_min_methods,
            exclude_methods=robust_exclude_methods,
            output_path=output_path,
        )
        summary_df = summarize_robust_sags(
            robust_df,
            output_path=output_path,
        )
        return {
            "per_category": marker_dict,
            "robust_sags": robust_df,
            "robust_sags_summary": summary_df,
        }

    return marker_dict


def _format_results(
    markers: pd.DataFrame,
    fmt: Literal["names", "z", "lfc", "df"],
) -> Union[GeneList, GeneZScores, GeneLFC, pd.DataFrame]:
    if fmt == "df":
        return markers
    elif fmt == "z":
        return list(zip(markers["names"], markers["z"]))
    elif fmt == "lfc":
        return list(zip(markers["names"], markers["logfoldchanges"]))
    else:
        return list(markers["names"].values)