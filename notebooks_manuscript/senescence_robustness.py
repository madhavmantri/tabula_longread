"""Robustness / control statistics for the senescence section (reviewer response).

Reusable pieces shared by the 06_36+ robustness notebooks. Everything here is
deliberately estimator-level and matrix-free where possible, so the cheap
obs-level checks can run in minutes without touching a 200k x 478k matrix.

Grouped by the reviewer point each function answers:

  #1 depth confound   depth_bin, depth_matched_indices, depth_stratified_test
  #2 donor as unit    per_donor_effects, dersimonian_laird, sign_consistency_p
  #5 power/equivalence  tost_equivalence, min_detectable_effect
  #6 co-detection     mantel_haenszel_or, depth_conditional_codetection,
                      matched_pair_or_null
  minor               cliffs_delta, bh, binom_direction_p

Nothing here reads project paths; callers pass arrays. That keeps it testable
and reusable by the Tier 2/3 notebooks.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

__all__ = [
    "bh", "cliffs_delta", "binom_direction_p",
    "depth_bin", "depth_bin_within", "depth_matched_indices",
    "matched_indices_by_block", "depth_stratified_test",
    "stratified_unit_effects", "rarefy_reference_mask",
    "per_donor_effects", "dersimonian_laird", "nested_meta_analysis",
    "cluster_t_summary", "sign_flip_test", "sign_consistency_p",
    "tost_equivalence", "min_detectable_effect",
    "mantel_haenszel_or", "depth_conditional_codetection", "matched_pair_or_null",
    "codetection_enrichment_by_unit",
    "read_h5ad_obs",
    # entropy machinery (must match 06_10b's metric exactly)
    "build_iso_to_gene", "counts_to_fractions", "per_cell_entropy",
    "per_cell_entropy_floor", "thin_to_depth", "chao_shen_entropy",
]


# ─────────────────────────────────────────────────────────────────────────────
# small shared utilities
# ─────────────────────────────────────────────────────────────────────────────
def bh(pvals):
    """Benjamini-Hochberg adjusted p-values. NaNs pass through as NaN."""
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    q = p[ok]
    n = q.size
    order = np.argsort(q)
    ranked = q[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    out[ok] = adj
    return out


def cliffs_delta(a, b):
    """Cliff's delta: P(a > b) - P(a < b), in [-1, 1].

    Rank-based so it pairs naturally with the Mann-Whitney tests already used
    in 06_10b. Computed from the MW U statistic in O(n log n) rather than the
    naive O(n*m) pairwise comparison, so it is usable at 200k cells.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    n, m = a.size, b.size
    if n == 0 or m == 0:
        return np.nan
    u = stats.mannwhitneyu(a, b, alternative="two-sided").statistic
    return float(2.0 * u / (n * m) - 1.0)


def binom_direction_p(n_higher, n_total, p_null=0.5):
    """Two-sided binomial p for 'effect went one way in n_higher of n_total'.

    The statistic the reviewer asks for on '19 of 21 testable cell types'.
    """
    if n_total == 0:
        return np.nan
    return float(stats.binomtest(int(n_higher), int(n_total), p_null,
                                 alternative="two-sided").pvalue)


def read_h5ad_obs(path, columns=None):
    """Read obs out of an .h5ad without loading X. Decodes categoricals.

    AnnData >=0.8 stores each obs column as its own dataset (or a group with
    categories/codes for categoricals); this reads them directly with h5py so a
    depth/donor/status table costs megabytes rather than the full object.
    """
    import h5py

    with h5py.File(path, "r") as f:
        g = f["obs"]
        order = g.attrs.get("column-order")
        if order is not None:
            cols = [c.decode() if isinstance(c, bytes) else c for c in order]
        else:
            cols = [k for k in g.keys() if not k.startswith("_")]
        idx_key = g.attrs.get("_index", "_index")
        if isinstance(idx_key, bytes):
            idx_key = idx_key.decode()

        def _read(name):
            node = g[name]
            if isinstance(node, h5py.Group):
                cats = [c.decode() if isinstance(c, bytes) else c
                        for c in node["categories"][:]]
                return pd.Categorical.from_codes(node["codes"][:], cats)
            arr = node[:]
            # h5py returns fixed-length strings as kind 'S' but VARIABLE-length
            # strings as dtype object holding bytes -- the index is usually the
            # latter, so decoding only on kind 'S' silently leaves b'...' keys
            # and every .reindex() against a str-indexed CSV joins to NaN.
            if arr.dtype.kind == "S":
                arr = np.array([x.decode() for x in arr])
            elif arr.dtype == object:
                arr = np.array([x.decode() if isinstance(x, bytes) else x
                                for x in arr])
            return arr

        want = cols if columns is None else [c for c in columns if c in cols]
        data = {c: _read(c) for c in want}
        index = _read(idx_key) if idx_key in g else None
    return pd.DataFrame(data, index=index)


# ─────────────────────────────────────────────────────────────────────────────
# #1  depth confound
# ─────────────────────────────────────────────────────────────────────────────
def depth_bin(depth, n_bins=10, log=True):
    """Assign cells to equal-frequency depth bins (quantile bins).

    Equal-frequency rather than equal-width because library size is heavily
    right-skewed; equal-width bins would put ~all cells in bin 0.
    """
    d = np.asarray(depth, dtype=float)
    x = np.log10(d + 1) if log else d
    ranks = stats.rankdata(x, method="average") / len(x)
    return np.clip((ranks * n_bins).astype(int), 0, n_bins - 1)


def depth_matched_indices(depth, group, n_bins=20, seed=0, ratio=1, within=None):
    """Depth-matched subsample: within each depth bin keep all minority-group
    cells and `ratio` x that many randomly drawn majority-group cells.

    Returns (idx_minority, idx_majority) into the original arrays. The two sets
    have, by construction, near-identical depth distributions -- so a difference
    that survives on them is not a depth artifact. This is the cheap version of
    the reviewer's ask (it matches WHICH cells are compared); the expensive
    version, downsampling the counts themselves, lives in the Tier 2 notebook.

    `within` optionally nests the match inside a grouping. Pass an array of
    labels (a tube id, a cell class, or a composite of both) and the whole
    procedure is run independently inside each label, with the depth bins
    recomputed as quantiles OF THAT GROUP rather than of the pooled cohort.
    Groups holding cells of only one arm contribute nothing.

    `within=None` (the default) reproduces the pooled behaviour exactly, so
    existing callers are unaffected. The distinction matters whenever the
    downstream estimator blocks on something: a pooled match balances depth
    only marginally, and a block that differs systematically in depth from the
    cohort can retain a within-block gap after it.

    NOTE on nesting: which arm is the minority is decided per group, not once
    globally. Pooled, the senescent arm is always the minority and every
    senescent cell survives; inside a thin group the reference arm can be
    scarcer, in which case reference cells are kept in full and senescent cells
    are subsampled to match. The 1:1 balance holds either way, but the
    "every senescent cell is retained" invariant does NOT, so callers that
    report senescent-cell counts should measure retention rather than assume
    it. Only the union of the two returned arrays is well defined as "the
    matched set"; do not assume the first return is the senescent arm.
    """
    if within is not None:
        w = np.asarray(within)
        g_all = np.asarray(group).astype(bool)
        d_all = np.asarray(depth, dtype=float)
        keep_min, keep_maj = [], []
        for lab in pd.unique(w):
            rows = np.where(w == lab)[0]
            if g_all[rows].sum() == 0 or (~g_all[rows]).sum() == 0:
                continue                      # nothing to match against inside this group
            a, c = depth_matched_indices(d_all[rows], g_all[rows], n_bins=n_bins,
                                         seed=seed, ratio=ratio, within=None)
            keep_min.append(rows[a])
            keep_maj.append(rows[c])
        if not keep_min:
            return np.array([], dtype=int), np.array([], dtype=int)
        return np.concatenate(keep_min), np.concatenate(keep_maj)

    rng = np.random.default_rng(seed)
    g = np.asarray(group).astype(bool)
    bins = depth_bin(depth, n_bins=n_bins)
    minority_is_true = g.sum() <= (~g).sum()
    keep_min, keep_maj = [], []
    for b in np.unique(bins):
        in_b = np.where(bins == b)[0]
        a = in_b[g[in_b]] if minority_is_true else in_b[~g[in_b]]
        c = in_b[~g[in_b]] if minority_is_true else in_b[g[in_b]]
        if a.size == 0 or c.size == 0:
            continue
        take = min(c.size, ratio * a.size)
        keep_min.append(a)
        keep_maj.append(rng.choice(c, take, replace=False))
    if not keep_min:
        return np.array([], dtype=int), np.array([], dtype=int)
    return np.concatenate(keep_min), np.concatenate(keep_maj)


def depth_bin_within(depth, unit, n_bins=10, log=True):
    """Per-cell depth bin computed as quantiles WITHIN each unit (tube/sample).

    `depth_bin` ranks the pooled cohort, so a cell from a shallow tube and a
    cell from a deep tube can share a global bin while sitting at opposite ends
    of their own tubes. Ranking inside the unit makes a bin mean "deep FOR THIS
    SAMPLE", which is what a within-sample depth control needs.

    Units are binned independently, so bin k in one unit is NOT the same depth
    range as bin k in another: always pair the returned label with the unit
    label when forming blocks (see `matched_indices_by_block`).
    """
    d = np.asarray(depth, dtype=float)
    u = np.asarray(unit).astype(str)
    out = np.zeros(d.size, dtype=int)
    for lab in pd.unique(u):
        rows = np.where(u == lab)[0]
        out[rows] = depth_bin(d[rows], n_bins=n_bins, log=log)
    return out


def matched_indices_by_block(group, block, seed=0, ratio=1):
    """BALANCED match inside each block. Returns (idx_group, idx_reference).

    Within each block both arms are capped, so exactly n_ref = ratio x n_group
    cells survive on each side and the returned sets are balanced BY
    CONSTRUCTION. Blocks holding only one arm contribute nothing.

    This is the difference from `depth_matched_indices`, which keeps the
    minority arm WHOLE (`take = min(c.size, ratio * a.size)`; the minority is
    appended unsubsampled). Pooled over a whole cohort the majority arm is
    almost always deep enough that the result comes out balanced anyway -- 06_45
    got exactly 3,134 vs 3,134. Inside a single tube it does not: a block can
    easily hold 4 senescent cells and 1 reference cell, and the unbalanced
    version would keep all 4. For an equal-n richness comparison (Fig. 4C) that
    silently breaks the very property the panel rests on, so callers doing
    within-sample matching should use this function instead.

    Unlike `depth_matched_indices`, the split is by `group`, not by whichever
    arm happens to be the minority, so the first return is always the True arm.

    Parameters
    ----------
    group : bool array, True for the arm to be preserved where possible
    block : label array; cells sharing a label are matched against each other.
            For within-sample depth matching pass unit + "|" + depth bin, with
            the bin from `depth_bin_within`.
    """
    rng = np.random.default_rng(seed)
    g = np.asarray(group).astype(bool)
    b = np.asarray(block).astype(str)
    r = max(1, int(ratio))
    order = np.arange(g.size)
    keep_g, keep_r = [], []
    for lab in np.unique(b):
        idx = order[b == lab]
        a = idx[g[idx]]
        c = idx[~g[idx]]
        if a.size == 0 or c.size == 0:
            continue
        n_a = min(a.size, c.size // r)
        if n_a == 0:                       # too few reference cells to pair even one
            continue
        keep_g.append(rng.choice(a, n_a, replace=False) if n_a < a.size else a)
        keep_r.append(rng.choice(c, r * n_a, replace=False))
    if not keep_g:
        return np.array([], dtype=int), np.array([], dtype=int)
    return np.concatenate(keep_g), np.concatenate(keep_r)


def depth_stratified_test(values, group, depth, n_bins=10):
    """Per-depth-bin Mann-Whitney + effect size, and a combined stratified test.

    Answers 'is the entropy difference constant across depth, or does it appear
    only in some depth range?' -- the supplementary panel the reviewer says
    would defuse point 1 entirely. The combined p uses van Elteren / stratified
    Wilcoxon weighting (sum of per-stratum z, weighted by stratum size).
    """
    v = np.asarray(values, dtype=float)
    g = np.asarray(group).astype(bool)
    bins = depth_bin(depth, n_bins=n_bins)
    rows, zs, ws = [], [], []
    for b in np.unique(bins):
        m = bins == b
        a, c = v[m & g], v[m & ~g]
        a, c = a[np.isfinite(a)], c[np.isfinite(c)]
        if a.size < 3 or c.size < 3:
            continue
        mw = stats.mannwhitneyu(a, c, alternative="two-sided")
        n1, n2 = a.size, c.size
        mu = n1 * n2 / 2.0
        sd = np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
        z = (mw.statistic - mu) / sd if sd > 0 else 0.0
        w = np.sqrt(n1 * n2 / (n1 + n2))
        rows.append({
            "depth_bin": int(b), "n_group": n1, "n_ref": n2,
            "median_group": float(np.median(a)), "median_ref": float(np.median(c)),
            "delta_median": float(np.median(a) - np.median(c)),
            "cliffs_delta": cliffs_delta(a, c), "p": float(mw.pvalue),
        })
        zs.append(z)
        ws.append(w)
    df = pd.DataFrame(rows)
    if zs:
        z_comb = np.dot(ws, zs) / np.sqrt(np.sum(np.square(ws)))
        p_comb = 2 * stats.norm.sf(abs(z_comb))
    else:
        z_comb, p_comb = np.nan, np.nan
    return df, float(z_comb), float(p_comb)


# ─────────────────────────────────────────────────────────────────────────────
# #2  donor as the unit of replication
# ─────────────────────────────────────────────────────────────────────────────
def per_donor_effects(values, group, donor, min_cells=20, unit_name="donor"):
    """Per-unit effect size + SE for a cell-level measurement.

    `unit_name` labels the output column, so the same routine serves donor-level
    and tube/sample-level replication (pass unit_name="tube_id"). NOTE the two
    are NOT interchangeable inferentially: tubes are nested within donors, so a
    meta-analysis over tubes treats correlated units as independent. Use
    `nested_meta_analysis` when the claim is about donors.

    Effect = difference in means (interpretable, and its SE is well defined);
    Cliff's delta reported alongside as the rank-based companion. Units with
    fewer than `min_cells` in either arm are returned with NaN so the caller can
    show them as untested rather than silently dropping them.
    """
    v = np.asarray(values, dtype=float)
    g = np.asarray(group).astype(bool)
    d = np.asarray(donor).astype(str)
    rows = []
    for dn in sorted(pd.unique(d)):
        m = d == dn
        a, c = v[m & g], v[m & ~g]
        a, c = a[np.isfinite(a)], c[np.isfinite(c)]
        row = {unit_name: dn, "n_group": a.size, "n_ref": c.size}
        if a.size >= min_cells and c.size >= min_cells:
            se = np.sqrt(a.var(ddof=1) / a.size + c.var(ddof=1) / c.size)
            row.update({
                "mean_group": a.mean(), "mean_ref": c.mean(),
                "effect": a.mean() - c.mean(), "se": se,
                "cliffs_delta": cliffs_delta(a, c),
                "p": stats.mannwhitneyu(a, c, alternative="two-sided").pvalue,
            })
        else:
            row.update({k: np.nan for k in
                        ("mean_group", "mean_ref", "effect", "se", "cliffs_delta", "p")})
        rows.append(row)
    return pd.DataFrame(rows)


def stratified_unit_effects(values, group, unit, strata, min_per_stratum=3,
                            unit_name="unit"):
    """Per-unit effect from a depth-STRATIFIED contrast that keeps EVERY
    comparison cell.

    1:1 depth matching controls the confound by throwing the comparison arm away
    -- on this dataset it discards ~148k of 151k non-senescent cells, and the
    surviving analysis then has to be thinned further, which injects noise and
    biases the contrast toward the null. Stratification achieves the same
    confound control without that loss: within each (unit, depth stratum) the
    senescent cells are contrasted against ALL reference cells at that depth,
    and the per-stratum differences are combined by inverse variance.

        effect_u = sum_s w_s d_s / sum_s w_s      w_s = 1 / (var_a/n_a + var_c/n_c)
        se_u     = sqrt(1 / sum_s w_s)

    A stratum contributes only if BOTH arms hold >= min_per_stratum cells, so
    strata where no comparison is possible drop out rather than being imputed.

    CAVEAT: this removes CONFOUNDING by depth (which cells are compared). It does
    not remove the depth-dependence of the entropy estimator itself within a
    stratum -- that is what Chao-Shen and count-thinning address. The two are
    complementary controls, not substitutes.
    """
    v = np.asarray(values, dtype=float)
    g = np.asarray(group).astype(bool)
    u = np.asarray(unit).astype(str)
    s = np.asarray(strata)
    rows = []
    for un in sorted(pd.unique(u)):
        m_u = u == un
        num_a = num_c = den = 0.0
        n_g = n_r = n_s = 0
        for st in np.unique(s[m_u]):
            m = m_u & (s == st)
            a, c = v[m & g], v[m & ~g]
            a, c = a[np.isfinite(a)], c[np.isfinite(c)]
            if a.size < min_per_stratum or c.size < min_per_stratum:
                continue
            var = a.var(ddof=1) / a.size + c.var(ddof=1) / c.size
            if not np.isfinite(var) or var <= 0:
                continue
            w = 1.0 / var
            num_a += w * a.mean()
            num_c += w * c.mean()
            den += w
            n_g += a.size
            n_r += c.size
            n_s += 1
        row = {unit_name: un, "n_strata": n_s, "n_group": n_g, "n_ref": n_r}
        if den > 0:
            # The two stratum-weighted arm means are returned as well as their
            # difference, so a paired per-unit plot can be drawn whose GAP IS
            # EXACTLY the effect used downstream:
            #   sum_s w_s (a_s - c_s) / sum_s w_s
            #     == sum_s w_s a_s / sum_s w_s  -  sum_s w_s c_s / sum_s w_s
            # Same weights, so the panel and the meta-analysis cannot disagree.
            # NOTE the weights are unit-specific, so these means are depth-
            # comparable WITHIN a unit; absolute level still differs between
            # units for real (tissue) reasons.
            row["mean_group_adj"] = float(num_a / den)
            row["mean_ref_adj"] = float(num_c / den)
            row["effect"] = float((num_a - num_c) / den)
            row["se"] = float(np.sqrt(1.0 / den))
        else:
            row["mean_group_adj"] = np.nan
            row["mean_ref_adj"] = np.nan
            row["effect"] = np.nan
            row["se"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def rarefy_reference_mask(group, unit, strata, seed=0, ratio=1):
    """Keep every group cell and `ratio` x that many reference cells, drawn
    WITHIN each (unit, stratum) block.

    Group-size control for stratified contrasts. Note this is a confirmation,
    not a correction: a difference of per-cell MEANS is unbiased under unequal
    arm sizes, since E[mean_a - mean_c] = mu_a - mu_c whatever n_a and n_c are.
    Unequal n changes precision only. That is unlike a POOLED richness measure
    (detected-isoform counts, aggregate per-gene entropy), which grows with n and
    genuinely does need rarefying -- see 06_10b Section 11, where rarefaction
    left the per-cell median stable while collapsing the aggregate-entropy gap.

    Rarefying within (unit, stratum) rather than globally preserves the blocking,
    so the rarefied estimate stays comparable to the unrarefied one.
    """
    rng = np.random.default_rng(seed)
    g = np.asarray(group).astype(bool)
    u = np.asarray(unit).astype(str)
    s = np.asarray(strata)
    keep = g.copy()
    key = np.array([a + "|" + str(b) for a, b in zip(u, s)])
    order = np.arange(g.size)
    for k in np.unique(key):
        idx = order[key == k]
        n_g = int(g[idx].sum())
        ref = idx[~g[idx]]
        if n_g == 0 or ref.size == 0:
            continue
        take = min(ref.size, ratio * n_g)
        keep[rng.choice(ref, take, replace=False)] = True
    return keep


def dersimonian_laird(effects, ses, knapp_hartung=True):
    """Random-effects meta-analysis (DerSimonian-Laird), Knapp-Hartung by default.

    Treats the supplied unit (donor, tube) as the unit of replication, so the
    pooled p cannot inherit the 150k-cell sample size that produces p = 5e-289.

    KNAPP-HARTUNG: plain DL computes SE = sqrt(1/sum w*) and tests against a
    NORMAL, which is only valid if tau^2 and the weights are KNOWN. They are
    estimated from the same handful of units, so the SE is too small and the
    test is anti-conservative -- with few units plus real heterogeneity a
    nominal 5% test can reject 10-20% of the time. Knapp-Hartung rescales the
    variance by the observed weighted dispersion,

        q     = 1/(k-1) * sum w*_i (e_i - pooled)^2
        se_hk = sqrt(q / sum w*)

    and refers it to t_{k-1}. `q` is floored at 1 so the correction can only
    ever WIDEN the interval -- without that guard, unusually homogeneous units
    drive q below 1 and produce an interval narrower than DL, the same
    pathology documented in `nested_meta_analysis`.

    With knapp_hartung=True (default) `effect`/`se`/`ci_*`/`p` are the corrected
    values and the uncorrected ones are returned under `z_*`. Set False to get
    the classic DL numbers as primary.
    """
    e = np.asarray(effects, dtype=float)
    s = np.asarray(ses, dtype=float)
    ok = np.isfinite(e) & np.isfinite(s) & (s > 0)
    e, s = e[ok], s[ok]
    k = e.size
    if k == 0:
        return {"k": 0}
    if k == 1:
        return {"k": 1, "effect": float(e[0]), "se": float(s[0]),
                "ci_low": float(e[0] - 1.96 * s[0]), "ci_high": float(e[0] + 1.96 * s[0]),
                "p": float(2 * stats.norm.sf(abs(e[0] / s[0]))),
                "Q": np.nan, "I2": np.nan, "tau2": 0.0, "method": "single unit"}
    w = 1.0 / s**2
    fixed = np.sum(w * e) / np.sum(w)
    Q = float(np.sum(w * (e - fixed) ** 2))
    df = k - 1
    C = np.sum(w) - np.sum(w**2) / np.sum(w)
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    w_star = 1.0 / (s**2 + tau2)
    pooled = float(np.sum(w_star * e) / np.sum(w_star))

    se_z = float(np.sqrt(1.0 / np.sum(w_star)))
    out = {
        "k": int(k), "Q": Q, "p_Q": float(stats.chi2.sf(Q, df)),
        "I2": float(max(0.0, (Q - df) / Q * 100) if Q > 0 else 0.0),
        "tau2": float(tau2),
        "z_effect": pooled, "z_se": se_z,
        "z_ci_low": pooled - 1.96 * se_z, "z_ci_high": pooled + 1.96 * se_z,
        "z_p": float(2 * stats.norm.sf(abs(pooled / se_z))),
    }

    q_raw = float(np.sum(w_star * (e - pooled) ** 2) / df)
    q_used = max(q_raw, 1.0)                       # correction may only widen
    se_hk = float(np.sqrt(q_used / np.sum(w_star)))
    tcrit = float(stats.t.isf(0.025, df))
    t_stat = pooled / se_hk if se_hk > 0 else np.inf
    out.update({
        "hk_q_raw": q_raw, "hk_q_used": q_used, "hk_se": se_hk, "hk_df": int(df),
        "hk_t": float(t_stat),
        "hk_ci_low": pooled - tcrit * se_hk, "hk_ci_high": pooled + tcrit * se_hk,
        "hk_p": float(2 * stats.t.sf(abs(t_stat), df)),
    })

    prefix = "hk_" if knapp_hartung else "z_"
    out["effect"] = pooled
    out["se"] = out[prefix + "se"]
    out["ci_low"] = out[prefix + "ci_low"]
    out["ci_high"] = out[prefix + "ci_high"]
    out["p"] = out[prefix + "p"]
    out["method"] = ("DerSimonian-Laird + Knapp-Hartung (t_{k-1})" if knapp_hartung
                     else "DerSimonian-Laird (z)")
    return out


def sign_flip_test(values, n_perm=200000, seed=0, exact_max_k=24):
    """Exact sign-flip permutation test that the values are centred on zero.

    Assumption-free alternative to the t-test on unit-level effects: under the
    null each unit's effect is equally likely to carry either sign, so the null
    distribution is generated by flipping signs. No normality assumption, which
    matters because at k=4 normality is both load-bearing and uncheckable.

    Exhaustive for k <= exact_max_k (all 2^k assignments), sampled above it.

    IMPORTANT FLOOR: the smallest attainable two-sided p is 2 / 2^k -- 0.125 at
    k=4, 4.8e-7 at k=22. Any p below that floor cannot come from an exact test
    on this many units; it can only come from a parametric assumption.
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    k = v.size
    floor = 2.0 ** (1 - k) if k > 0 else np.nan
    if k == 0:
        return {"k": 0, "p": np.nan, "floor": np.nan}
    obs = abs(float(v.mean()))

    if k <= exact_max_k:
        # CHUNKED enumeration. Building all 2^k x k sign vectors at once needs
        # ~740 MB of int64 intermediates at k=22, which is why an earlier
        # version capped exact_max_k at 16 and silently fell through to
        # sampling -- returning p = 1/(n_perm+1), the resolution limit, rather
        # than a measurement. Chunking makes k up to ~24 exact and cheap.
        total = 1 << k
        bits = np.arange(k, dtype=np.int64)
        hits = 0
        chunk = 1 << 18
        for start in range(0, total, chunk):
            idx = np.arange(start, min(start + chunk, total), dtype=np.int64)
            signs = 1 - 2 * ((idx[:, None] >> bits) & 1).astype(np.int8)
            means = (signs * v).mean(axis=1)
            hits += int((np.abs(means) >= obs - 1e-12).sum())
        p = hits / total
        exact = True
        n_used = total
        resolution = 1.0 / total
    else:
        rng = np.random.default_rng(seed)
        hits = 0
        done = 0
        chunk = 20000
        while done < n_perm:
            m = min(chunk, n_perm - done)
            signs = rng.choice(np.array([-1.0, 1.0]), size=(m, k))
            means = (signs * v).mean(axis=1)
            hits += int((np.abs(means) >= obs - 1e-12).sum())
            done += m
        p = float(hits + 1) / (n_perm + 1)
        exact = False
        n_used = n_perm
        resolution = 1.0 / (n_perm + 1)
    return {"k": int(k), "observed_mean": float(v.mean()), "p": p,
            "exact": exact, "n_assignments": int(n_used), "floor": float(floor),
            "resolution": float(resolution),
            "resolution_limited": bool(not exact and hits == 0),
            "n_positive": int((v > 0).sum())}


def cluster_t_summary(values):
    """One-sample t summary treating each supplied value as ONE observation.

    Use when the k cluster-level estimates ARE the sample -- e.g. four donor
    effects standing in for four people. Deliberately ignores the precision of
    each individual estimate, because that precision is driven by cell counts
    and would otherwise smuggle cell-level replication back into what is meant
    to be a person-level claim.
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    k = v.size
    if k < 2:
        return {"k": int(k), "effect": float(v[0]) if k else np.nan,
                "sd": np.nan, "se": np.nan, "ci_low": np.nan, "ci_high": np.nan,
                "t": np.nan, "df": max(k - 1, 0), "p": np.nan}
    mean = float(v.mean())
    sd = float(v.std(ddof=1))
    se = sd / np.sqrt(k)
    df = k - 1
    tcrit = float(stats.t.isf(0.025, df))
    t = mean / se if se > 0 else np.inf
    return {"k": int(k), "effect": mean, "sd": sd, "se": float(se),
            "ci_low": mean - tcrit * se, "ci_high": mean + tcrit * se,
            "t": float(t), "df": int(df),
            "p": float(2 * stats.t.sf(abs(t), df)) if se > 0 else 0.0}


def nested_meta_analysis(unit_effects, unit_ses, unit_parent):
    """Two-stage analysis for units nested inside a higher-level cluster.

    Tubes are nested within donors: several tubes from one person share genotype,
    tissue handling and library prep, so pooling tubes as if independent
    understates uncertainty for any donor-level claim.

      stage 1  pool units WITHIN each parent (random effects) -> one estimate
               per parent
      stage 2  summarise the parent estimates as k observations (one-sample t)

    WHY STAGE 2 IS A t-TEST AND NOT ANOTHER DerSimonian-Laird POOL: with few and
    mutually concordant parents DL estimates tau^2 ~ 0, the random-effects model
    collapses to fixed effects, and the pooled SE becomes 1/sqrt(sum 1/se_i^2) --
    which inherits the CELL-driven precision of stage 1 and reintroduces exactly
    the pseudoreplication the nesting was meant to remove. Measured on this
    dataset: a DL stage 2 over k=4 donors returned p = 1.8e-24 with a CI NARROWER
    than the flat over-tubes pooling, which is diagnostic of the failure. The t
    summary over the same four donor estimates gives p ~ 2e-3.

    DL numbers are still returned under `dl_*` keys for comparison, but
    `effect`/`ci_low`/`ci_high`/`p` are the t-based ones -- quote those.
    Returns (stage2_summary, per_parent_dataframe).
    """
    e = np.asarray(unit_effects, dtype=float)
    s = np.asarray(unit_ses, dtype=float)
    p = np.asarray(unit_parent).astype(str)
    ok = np.isfinite(e) & np.isfinite(s) & (s > 0)
    e, s, p = e[ok], s[ok], p[ok]

    rows = []
    for parent in sorted(pd.unique(p)):
        m = p == parent
        sub = dersimonian_laird(e[m], s[m])
        rows.append({"parent": parent, "n_units": int(m.sum()),
                     "effect": sub.get("effect", np.nan),
                     "se": sub.get("se", np.nan),
                     "I2_within": sub.get("I2", np.nan)})
    per_parent = pd.DataFrame(rows)

    stage2 = cluster_t_summary(per_parent["effect"].values)
    # assumption-free companion: at small k the t p-value is carried entirely by
    # a normality assumption that cannot be checked, so report the exact
    # sign-flip p and its floor next to it.
    sf = sign_flip_test(per_parent["effect"].values)
    stage2["signflip_p"] = sf["p"]
    stage2["signflip_floor"] = sf["floor"]
    stage2["signflip_exact"] = sf["exact"]
    # knapp_hartung=False ON PURPOSE: the dl_* keys exist to DEMONSTRATE the
    # pathology (z-test + tau^2 -> 0 inheriting stage-1 cell-level precision).
    # With the KH correction on, that failure no longer occurs, so a corrected
    # comparison would silently stop showing the thing it is there to show.
    dl = dersimonian_laird(per_parent["effect"].values, per_parent["se"].values,
                           knapp_hartung=False)
    for key in ("effect", "se", "ci_low", "ci_high", "p", "I2", "tau2"):
        if key in dl:
            stage2["dl_" + key] = dl[key]
    # ...and the KH-corrected version alongside, which SHOULD behave sensibly
    dl_hk = dersimonian_laird(per_parent["effect"].values, per_parent["se"].values,
                              knapp_hartung=True)
    stage2["dl_hk_p"] = dl_hk.get("p", np.nan)
    stage2["dl_hk_ci_low"] = dl_hk.get("ci_low", np.nan)
    stage2["dl_hk_ci_high"] = dl_hk.get("ci_high", np.nan)
    stage2["stage2_method"] = "one-sample t over parent estimates"
    return stage2, per_parent


def sign_consistency_p(effects):
    """Binomial p that all donors agree in sign (the '4/4 donors' statistic).

    The reviewer's point: consistency across the 4 biological replicates is more
    persuasive than a p-value with 291 zeros. With k=4 the floor is p=0.125,
    which is worth stating plainly rather than dressing up.
    """
    e = np.asarray(effects, dtype=float)
    e = e[np.isfinite(e)]
    if e.size == 0:
        return {"k": 0, "n_positive": 0, "p": np.nan}
    n_pos = int((e > 0).sum())
    n_high = max(n_pos, e.size - n_pos)
    return {"k": int(e.size), "n_positive": n_pos,
            "all_same_sign": bool(n_pos in (0, e.size)),
            "p": binom_direction_p(n_high, e.size),
            "p_floor_note": f"minimum attainable two-sided p at k={e.size} is "
                            f"{binom_direction_p(e.size, e.size):.3g}"}


# ─────────────────────────────────────────────────────────────────────────────
# #5  power / equivalence
# ─────────────────────────────────────────────────────────────────────────────
def tost_equivalence(a, b, margin):
    """Two one-sided tests: is |mean(a) - mean(b)| smaller than `margin`?

    Converts 'not significant' into 'significantly null', which is what the
    reviewer wants before a cell class is called unaffected. p < 0.05 means
    statistically equivalent within +/- margin.
    """
    a = np.asarray(a, dtype=float); a = a[np.isfinite(a)]
    b = np.asarray(b, dtype=float); b = b[np.isfinite(b)]
    if a.size < 3 or b.size < 3:
        return {"n_a": a.size, "n_b": b.size, "p_tost": np.nan, "equivalent": None}
    diff = a.mean() - b.mean()
    se = np.sqrt(a.var(ddof=1) / a.size + b.var(ddof=1) / b.size)
    if se == 0:
        return {"n_a": a.size, "n_b": b.size, "p_tost": np.nan, "equivalent": None}
    dof = (se**2) ** 2 / (
        (a.var(ddof=1) / a.size) ** 2 / (a.size - 1)
        + (b.var(ddof=1) / b.size) ** 2 / (b.size - 1))
    p_lower = stats.t.sf((diff + margin) / se, dof)      # H0: diff <= -margin
    p_upper = stats.t.cdf((diff - margin) / se, dof)     # H0: diff >= +margin
    p_tost = float(max(p_lower, p_upper))
    return {"n_a": a.size, "n_b": b.size, "diff": float(diff), "se": float(se),
            "margin": float(margin), "p_tost": p_tost, "equivalent": p_tost < 0.05}


def min_detectable_effect(n_a, n_b, sd, alpha=0.05, power=0.8):
    """Smallest difference in means detectable at the given n -- the number that
    distinguishes 'cell-type-restricted' from 'underpowered'."""
    if n_a < 2 or n_b < 2 or not np.isfinite(sd) or sd <= 0:
        return np.nan
    z_a = stats.norm.isf(alpha / 2)
    z_b = stats.norm.isf(1 - power)
    return float((z_a + z_b) * sd * np.sqrt(1.0 / n_a + 1.0 / n_b))


# ─────────────────────────────────────────────────────────────────────────────
# #6  co-detection null
# ─────────────────────────────────────────────────────────────────────────────
def mantel_haenszel_or(x, y, strata):
    """Mantel-Haenszel odds ratio for x vs y, pooled over strata.

    Stratifying on capture depth is the direct answer to 'the OR is inflated by
    cell-level capture variation': if the association is purely a depth effect
    the MH OR collapses toward 1 while the crude OR stays large.
    """
    x = np.asarray(x).astype(bool)
    y = np.asarray(y).astype(bool)
    s = np.asarray(strata)
    num = den = 0.0
    rows = []
    for lev in np.unique(s):
        m = s == lev
        a = int((x[m] & y[m]).sum())
        b = int((x[m] & ~y[m]).sum())
        c = int((~x[m] & y[m]).sum())
        d = int((~x[m] & ~y[m]).sum())
        n = a + b + c + d
        if n == 0:
            continue
        num += a * d / n
        den += b * c / n
        rows.append({"stratum": lev, "a": a, "b": b, "c": c, "d": d,
                     "or": (a * d) / (b * c) if b * c > 0 else np.inf})
    mh = num / den if den > 0 else np.inf
    a_ = int((x & y).sum()); b_ = int((x & ~y).sum())
    c_ = int((~x & y).sum()); d_ = int((~x & ~y).sum())
    crude = (a_ * d_) / (b_ * c_) if b_ * c_ > 0 else np.inf
    return {"mh_or": float(mh), "crude_or": float(crude),
            "per_stratum": pd.DataFrame(rows)}


def depth_conditional_codetection(x, y, depth, n_bins=20, n_perm=1000, seed=0):
    """Observed co-detection vs a depth-conditional independence null.

    Null: within a depth stratum, p16 and p14ARF detection are independent.
    Implemented two ways, which should agree:
      - analytic   E[co] = sum over strata of n_s * p_s * q_s
      - permutation  shuffle y within each stratum, n_perm times
    A large excess over this null is co-regulation; no excess means the crude OR
    was capture-depth variation, exactly as the reviewer suspects.
    """
    rng = np.random.default_rng(seed)
    x = np.asarray(x).astype(bool)
    y = np.asarray(y).astype(bool)
    bins = depth_bin(depth, n_bins=n_bins)
    observed = int((x & y).sum())

    expected = 0.0
    for b in np.unique(bins):
        m = bins == b
        n = int(m.sum())
        if n == 0:
            continue
        expected += n * x[m].mean() * y[m].mean()

    null = np.empty(n_perm)
    for i in range(n_perm):
        yp = y.copy()
        for b in np.unique(bins):
            idx = np.where(bins == b)[0]
            yp[idx] = rng.permutation(y[idx])
        null[i] = (x & yp).sum()
    p_emp = float((null >= observed).sum() + 1) / (n_perm + 1)
    return {
        "observed": observed,
        "expected_analytic": float(expected),
        "expected_perm_mean": float(null.mean()),
        "perm_sd": float(null.std()),
        "enrichment": float(observed / expected) if expected > 0 else np.inf,
        "z": float((observed - null.mean()) / null.std()) if null.std() > 0 else np.inf,
        "p_perm": p_emp,
        "null": null,
    }


def codetection_enrichment_by_unit(x, y, unit, strata, unit_name="unit",
                                   min_observed=3):
    """Per-unit co-detection enrichment over a depth-conditional null.

    The sample-level companion to `depth_conditional_codetection`. For each unit
    (tube, donor) and each stratum within it, the expected co-detection under
    within-stratum independence is n_s * p(x)_s * p(y)_s; observed and expected
    are summed over strata and reported as log(O/E) with a standard log-ratio
    standard error.

    WHY log(O/E) AND NOT THE ODDS RATIO: at sample level the 2x2 is degenerate.
    On this dataset 20 of 28 tubes contain ZERO p14ARF-only cells, so the odds
    ratio is infinite or undefined in most samples and could only be
    meta-analysed after Haldane-Anscombe correcting the majority of strata --
    fitting noise. O/E needs only the co-detection count and the two marginals,
    so an empty off-diagonal cell does not break it.

    `min_observed` guards log(0); units below it return NaN so the caller shows
    them as untested rather than silently dropping them.
    """
    x = np.asarray(x).astype(bool)
    y = np.asarray(y).astype(bool)
    u = np.asarray(unit).astype(str)
    s = np.asarray(strata)
    rows = []
    for un in sorted(pd.unique(u)):
        m_u = u == un
        obs_n = 0
        exp_n = 0.0
        n_s = 0
        for st in np.unique(s[m_u]):
            m = m_u & (s == st)
            n = int(m.sum())
            if n == 0:
                continue
            px, py = x[m].mean(), y[m].mean()
            if px == 0 or py == 0:
                continue
            exp_n += n * px * py
            obs_n += int((x[m] & y[m]).sum())
            n_s += 1
        row = {unit_name: un, "n_strata": n_s, "n_cells": int(m_u.sum()),
               "observed": obs_n, "expected": float(exp_n)}
        if obs_n >= min_observed and exp_n > 0:
            row["enrichment"] = float(obs_n / exp_n)
            row["effect"] = float(np.log(obs_n / exp_n))       # log scale for pooling
            row["se"] = float(np.sqrt(1.0 / obs_n + 1.0 / exp_n))
        else:
            row["enrichment"] = np.nan
            row["effect"] = np.nan
            row["se"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def matched_pair_or_null(detect_matrix, target_rates, n_pairs=500, tol=0.1, seed=0):
    """Empirical null for the co-detection OR from unlinked gene pairs matched
    on marginal detection rate.

    `detect_matrix` is cells x genes boolean (or 0/1). For each draw, pick two
    genes whose detection rates are within `tol` (relative) of the two target
    rates and compute their OR. The resulting distribution says how large an OR
    a pair of *independent* genes at this detection rate produces purely from
    cell-level capture variation -- the reviewer's 'matched for expression'
    control. Genes are drawn from different chromosomes-worth of the matrix only
    in the sense that they are random, so callers should exclude the CDKN2A
    locus itself before passing the matrix in.
    """
    rng = np.random.default_rng(seed)
    D = detect_matrix
    rates = np.asarray(D.mean(axis=0)).ravel()
    r1, r2 = target_rates

    def _candidates(r):
        lo, hi = r * (1 - tol), r * (1 + tol)
        c = np.where((rates >= lo) & (rates <= hi))[0]
        return c

    c1, c2 = _candidates(r1), _candidates(r2)
    if c1.size == 0 or c2.size == 0:
        return {"n_pairs": 0, "note": "no genes matched the target detection rates"}
    ors = []
    for _ in range(n_pairs):
        i = rng.choice(c1)
        j = rng.choice(c2)
        if i == j:
            continue
        xi = np.asarray(D[:, i]).ravel().astype(bool)
        yj = np.asarray(D[:, j]).ravel().astype(bool)
        a = int((xi & yj).sum()); b = int((xi & ~yj).sum())
        c = int((~xi & yj).sum()); d = int((~xi & ~yj).sum())
        if b * c == 0:
            continue
        ors.append((a * d) / (b * c))
    ors = np.asarray(ors, dtype=float)
    if ors.size == 0:
        return {"n_pairs": 0, "note": "no evaluable pairs"}
    return {
        "n_pairs": int(ors.size),
        "median_or": float(np.median(ors)),
        "p95_or": float(np.percentile(ors, 95)),
        "p99_or": float(np.percentile(ors, 99)),
        "max_or": float(ors.max()),
        "ors": ors,
    }


# ─────────────────────────────────────────────────────────────────────────────
# entropy machinery
#
# These reimplement `_per_cell_entropy_batched` from 06_10b. The definitions MUST
# stay identical to it -- the whole point of the depth controls is to compare a
# corrected estimate against the published one, and that comparison is meaningless
# if the metric drifted. Verified against the cached values: the recomputed
# depth-matched delta (0.0362) matches the cache-derived one (0.0374) to within
# resampling noise.
# ─────────────────────────────────────────────────────────────────────────────
def build_iso_to_gene(gene_labels):
    """Isoform -> gene incidence matrix, plus the per-isoform gene index.

    Returns (iso_to_gene [n_iso x n_genes csc], gene_of_iso, unique_genes,
    gene_col_indices) -- the four objects every entropy routine below needs.
    """
    from scipy.sparse import csc_matrix

    genes = np.asarray(gene_labels).astype(str)
    unique_genes, gene_of_iso = np.unique(genes, return_inverse=True)
    iso_to_gene = csc_matrix(
        (np.ones(len(genes), np.float32), (np.arange(len(genes)), gene_of_iso)),
        shape=(len(genes), len(unique_genes)))
    gene_col_indices = [np.where(gene_of_iso == g)[0]
                        for g in range(len(unique_genes))]
    return iso_to_gene, gene_of_iso, unique_genes, gene_col_indices


def counts_to_fractions(X, iso_to_gene, gene_of_iso):
    """Per-cell isoform fractions within each gene, from raw counts."""
    from scipy.sparse import csr_matrix, issparse

    X = X.tocsr() if issparse(X) else csr_matrix(X)
    gene_tot = np.asarray((X @ iso_to_gene).todense())
    Xc = X.tocoo()
    denom = gene_tot[Xc.row, gene_of_iso[Xc.col]]
    data = np.divide(Xc.data, denom,
                     out=np.zeros(Xc.data.shape, dtype=np.float64),
                     where=denom > 0)
    return csr_matrix((data, (Xc.row, Xc.col)), shape=X.shape)


def per_cell_entropy(frac, iso_to_gene, min_isoforms_per_gene=2, batch_size=2000):
    """Mean over genes of the Shannon entropy of that gene's isoform fractions,
    counting only genes with >= min_isoforms_per_gene EXPRESSED isoforms.

    Identical to `_per_cell_entropy_batched` in 06_10b.
    """
    from scipy.sparse import issparse

    n_cells = frac.shape[0]
    out = np.zeros(n_cells, dtype=np.float64)
    for start in range(0, n_cells, batch_size):
        end = min(start + batch_size, n_cells)
        batch = frac[start:end]
        if issparse(batch):
            b = batch.tocsr()
            term = b.copy()
            term.data = -term.data * np.log2(np.maximum(term.data, 1e-300))
            gene_entropy = np.asarray((term @ iso_to_gene).todense())
            n_expr = np.asarray(((b > 0).astype(np.float32) @ iso_to_gene).todense())
        else:
            arr = np.asarray(batch)
            term = np.where(arr > 0, -arr * np.log2(np.maximum(arr, 1e-300)), 0.0)
            gene_entropy = np.asarray(term @ iso_to_gene)
            n_expr = np.asarray((arr > 0).astype(np.float32) @ iso_to_gene)
        valid = n_expr >= min_isoforms_per_gene
        s = np.where(valid, gene_entropy, 0.0).sum(axis=1)
        k = valid.sum(axis=1)
        out[start:end] = np.divide(s, k, out=np.zeros(s.shape, dtype=np.float64),
                                   where=k > 0)
    return out


def per_cell_entropy_floor(frac, gene_tot, iso_to_gene, min_counts,
                           min_isoforms_per_gene=2, batch_size=2000):
    """Per-cell entropy counting only genes with >= min_counts IN THAT CELL.

    The reviewer asks to "restrict to genes with >= N counts in every cell". A
    gene-level floor across all cells is unsatisfiable in single-cell data (0 of
    26,940 genes reach 5 counts in every one of ~6k cells), so the floor is
    applied per (cell, gene): no entropy value is ever estimated from 2 UMIs.
    Returns (entropy, n_genes_used_per_cell).
    """
    n_cells = frac.shape[0]
    out = np.full(n_cells, np.nan)
    n_used = np.zeros(n_cells, dtype=int)
    for start in range(0, n_cells, batch_size):
        end = min(start + batch_size, n_cells)
        b = frac[start:end].tocsr()
        term = b.copy()
        term.data = -term.data * np.log2(np.maximum(term.data, 1e-300))
        gene_entropy = np.asarray((term @ iso_to_gene).todense())
        n_expr = np.asarray(((b > 0).astype(np.float32) @ iso_to_gene).todense())
        valid = (n_expr >= min_isoforms_per_gene) & (gene_tot[start:end] >= min_counts)
        s = np.where(valid, gene_entropy, 0.0).sum(axis=1)
        k = valid.sum(axis=1)
        out[start:end] = np.divide(s, k, out=np.full(s.shape, np.nan), where=k > 0)
        n_used[start:end] = k
    return out, n_used


def thin_to_depth(X, target, rng):
    """Binomial thinning of each cell's counts down to `target` total.

    p_i = min(1, target / rowsum_i); each count c -> Binomial(c, p_i). Preserves
    the multinomial structure, which is what makes the downsampled entropy
    directly comparable across cells of originally different depth.
    """
    from scipy.sparse import csr_matrix

    X = X.tocsr()
    rowsum = np.asarray(X.sum(axis=1)).ravel()
    p = np.where(rowsum > 0, np.minimum(1.0, target / np.maximum(rowsum, 1)), 0.0)
    Xc = X.tocoo()
    newdata = rng.binomial(Xc.data.astype(np.int64), p[Xc.row])
    keep = newdata > 0
    return csr_matrix((newdata[keep], (Xc.row[keep], Xc.col[keep])), shape=X.shape)


def chao_shen_entropy(X, gene_col_indices, min_isoforms=2):
    """Chao-Shen bias-corrected per-cell entropy (mean over qualifying genes).

    Plug-in Shannon UNDERestimates entropy at small counts, and the size of that
    bias depends on per-cell depth -- which is why the plug-in estimate alone
    cannot separate biology from sequencing depth. Chao-Shen applies a
    Good-Turing coverage correction:
        C = 1 - f1/N,  p_a = C * p_hat,  H = -sum p_a log p_a / (1 - (1-p_a)^N)
    Returns (entropy, n_genes_used_per_cell).
    """
    X = X.tocsc()
    n_cells = X.shape[0]
    tot = np.zeros(n_cells)
    cnt = np.zeros(n_cells, dtype=int)
    for cols in gene_col_indices:
        if len(cols) < min_isoforms:
            continue
        sub = X[:, cols].toarray()
        N = sub.sum(axis=1)
        nz = (sub > 0).sum(axis=1)
        ok = (N > 0) & (nz >= min_isoforms)
        if not ok.any():
            continue
        s = sub[ok]
        Nn = N[ok][:, None]
        f1 = (s == 1).sum(axis=1)[:, None]
        C = 1.0 - f1 / Nn
        C = np.where(C <= 0, 1.0 / Nn, C)          # guard all-singleton cells
        pa = C * (s / Nn)
        with np.errstate(divide="ignore", invalid="ignore"):
            denom = 1.0 - np.power(1.0 - pa, Nn)
            term = np.where(pa > 0,
                            -pa * np.log2(pa) / np.where(denom > 0, denom, 1), 0.0)
        idx = np.where(ok)[0]
        tot[idx] += term.sum(axis=1)
        cnt[idx] += 1
    return np.divide(tot, cnt, out=np.zeros(tot.shape), where=cnt > 0), cnt
