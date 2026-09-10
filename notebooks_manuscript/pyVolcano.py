import inspect
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text


def _adjust_text_safe(texts, **kwargs):
    """Call adjust_text with only the kwargs this installed version accepts.

    adjustText renamed/removed several arguments between 0.8 and 1.x (expand_text ->
    expand, force_points -> force_static, ...). Filtering against the live signature
    keeps this helper working across versions instead of raising TypeError.
    """
    params = inspect.signature(adjust_text).parameters
    accepts_var_kw = any(p.kind is inspect.Parameter.VAR_KEYWORD for p in params.values())
    kw = {k: v for k, v in kwargs.items() if v is not None and (accepts_var_kw or k in params)}
    dropped = [k for k, v in kwargs.items()
               if v is not None and not (accepts_var_kw or k in params)]
    if dropped:
        print(f"[pyVolcano] adjust_text ignored unsupported kwargs: {dropped}")
    return adjust_text(texts, **kw)

def plot_volcano(
    dge_results_df, 
    lfc_threshold = 0.5, 
    pval_threshold = 0.01,
    group=None, 
    pval_col ='pvals_adj', 
    logfc_col ='logfoldchanges', 
    name_col ='names',
    title='Volcano Plot',
    xlabel=None,
    all_text_size = None,
    text_size = 6,
    title_text_size = 6,
    legend_text_size = 6,
    legend_loc = 6,
    legend_bbox_to_anchor = (1.0, 0.5),
    axis_title_text_size = 6,
    axis_tick_text_size = 6,
    xlims=None, 
    ylims=None,
    savefig=None, 
    show_grid = False,
    point_alpha = 0.5,
    point_size = 10, 
    figsize=(10, 6),
    
    # --- NEW PARAMETERS -----------------------------------------------------
    show_gene_labels   = True,        # NEW: toggle text labels
    label_lfc_threshold= None,        # NEW: lfc threshold for labels
    label_pval_threshold=None,        # NEW: p-value threshold for labels
    return_gene_lists  = True,        # NEW: toggle return of gene lists
    show_legend        = True,        # NEW: hide significance legend
    # ------------------------------------------------------------------------

    # --- label de-crowding (passed through to adjustText) -------------------
    # Defaults reproduce adjustText's own behaviour, so existing figures are
    # unchanged unless a caller opts in. The single most effective knob is
    # max_labels; expand/force_text spread out whatever labels remain.
    label_genes        = None,        # explicit gene names to label, thresholds ignored
    label_genes_only   = True,        # True: label ONLY label_genes; False: add to the
                                      #       threshold-selected set
    max_labels         = None,        # keep only the N most significant labels
    label_expand       = None,        # e.g. (1.6, 2.0) — bbox padding per label
    label_force_text   = None,        # e.g. (0.5, 1.0) — label<->label repulsion
    label_force_static = None,        # e.g. (0.3, 0.5) — label<->point repulsion
    label_min_arrow_len= None,        # suppress leader lines shorter than this
    adjust_text_kwargs = None,        # escape hatch; merged last, wins ties
    # ------------------------------------------------------------------------
):
    """
    Plot a volcano plot for differential expression analysis results.

    Returns
    -------
    (down_list, up_list) if return_gene_lists is True, else None.
    """
    # ---------- defaults for the new thresholds -----------------------------
    if label_lfc_threshold is None:
        label_lfc_threshold = lfc_threshold
    if label_pval_threshold is None:
        label_pval_threshold = pval_threshold
    # ------------------------------------------------------------------------

    if all_text_size is not None:
        text_size = title_text_size = legend_text_size = \
                    axis_title_text_size = axis_tick_text_size = all_text_size

    df              = dge_results_df
    pvals           = df[pval_col]
    logfoldchanges  = df[logfc_col]
    gene_names      = df[name_col]

    sig_pval, sig_logfc = pval_threshold, lfc_threshold

    plt.figure(figsize=figsize)
    plt.grid(show_grid)
    plt.rcParams.update({'font.size': text_size})

    # base layer
    plt.scatter(logfoldchanges, -np.log10(pvals),
                color='grey', alpha=point_alpha, s=point_size, label='NS')

    # significance layers
    sig_pts = (pvals < sig_pval) & (np.abs(logfoldchanges) > sig_logfc)
    plt.scatter(logfoldchanges[sig_pts], -np.log10(pvals[sig_pts]),
                color='red', alpha=point_alpha, s=point_size,
                label='p-value & log2 FC')

    pv_only = (pvals < sig_pval) & ~sig_pts
    plt.scatter(logfoldchanges[pv_only], -np.log10(pvals[pv_only]),
                color='blue', alpha=point_alpha, s=point_size, label='p-value')

    lfc_only = (~(pvals < sig_pval)) & (np.abs(logfoldchanges) > sig_logfc)
    plt.scatter(logfoldchanges[lfc_only], -np.log10(pvals[lfc_only]),
                color='green', alpha=point_alpha, s=point_size, label='log2 FC')

    # reference lines
    plt.axhline(-np.log10(sig_pval), color='black', linestyle='--')
    plt.axvline( sig_logfc,  color='black', linestyle='--')
    plt.axvline(-sig_logfc,  color='black', linestyle='--')
    plt.tick_params(axis='both', which='major', pad=0.2, size=5)

    # ---------------- label logic ------------------------------------------
    texts, down_reg, up_reg = [], [], []
    label_mask = (pvals < label_pval_threshold) & \
                 (np.abs(logfoldchanges) > label_lfc_threshold)

    label_idx = np.where(label_mask)[0]

    # De-crowding: keep only the most significant labels. The gene LISTS returned to
    # the caller are built from the full label_mask below, so capping affects what is
    # DRAWN, never what is reported.
    draw_idx = label_idx
    if max_labels is not None and len(label_idx) > max_labels:
        _p = np.asarray(pvals)[label_idx]
        _l = np.asarray(logfoldchanges)[label_idx]
        score = -np.log10(np.maximum(_p, 1e-300)) * np.abs(_l)
        # split the budget across both arms so a one-sided volcano is not produced
        # by the capping itself
        up_pos = label_idx[_l > 0][np.argsort(-score[_l > 0])]
        dn_pos = label_idx[_l < 0][np.argsort(-score[_l < 0])]
        half = int(np.ceil(max_labels / 2))
        keep = list(up_pos[:half]) + list(dn_pos[:half])
        if len(keep) > max_labels:                       # trim the weaker arm
            keep = list(np.array(keep)[np.argsort(
                -np.array([-np.log10(max(pvals[k], 1e-300)) * abs(logfoldchanges[k])
                           for k in keep]))][:max_labels])
        draw_idx = np.array(sorted(keep), dtype=int)
        print(f"[pyVolcano] labelling {len(draw_idx)} of {len(label_idx)} "
              f"significant points (max_labels={max_labels})")

    # Explicit gene list: label exactly these, regardless of any threshold. Applied
    # AFTER the max_labels cap so a hand-picked gene is never silently dropped.
    if label_genes is not None:
        want = {str(g) for g in label_genes}
        gname = np.asarray(gene_names).astype(str)
        explicit_idx = np.where(np.isin(gname, list(want)))[0]
        missing = sorted(want - set(gname[explicit_idx]))
        if missing:
            print(f"[pyVolcano] label_genes not found in '{name_col}': {missing}")
        draw_idx = (explicit_idx if label_genes_only
                    else np.union1d(np.asarray(draw_idx, dtype=int), explicit_idx))
        print(f"[pyVolcano] labelling {len(draw_idx)} gene(s) from label_genes "
              f"({'only these' if label_genes_only else 'plus the threshold set'})")

    # Gene LISTS returned to the caller stay threshold-based and are unaffected by
    # max_labels / label_genes, which control only what is DRAWN.
    if return_gene_lists:
        for i in label_idx:
            if logfoldchanges[i] < -label_lfc_threshold:
                down_reg.append((gene_names[i], logfoldchanges[i]))
            elif logfoldchanges[i] >  label_lfc_threshold:
                up_reg.append((gene_names[i], logfoldchanges[i]))

    if show_gene_labels:                                            # CHANGED
        for i in draw_idx:
            texts.append(
                plt.text(logfoldchanges[i], -np.log10(pvals[i]),
                         gene_names[i], fontsize=text_size,
                         ha='center', va='center'))

    if show_gene_labels and texts:                                   # CHANGED
        _kw = dict(arrowprops=dict(arrowstyle='-', color='black'),
                   expand=label_expand,
                   force_text=label_force_text,
                   force_static=label_force_static,
                   min_arrow_len=label_min_arrow_len)
        if adjust_text_kwargs:
            _kw.update(adjust_text_kwargs)
        _adjust_text_safe(texts, **_kw)
    # -----------------------------------------------------------------------

    plt.xlabel(xlabel if xlabel is not None else r'$log_2$ fold-change',
               fontsize=axis_title_text_size, labelpad=0.1)

    # CHANGED: mathtext-safe ylabel (no \text command)
    # BOLD: fontweight alone cannot bold a string that is entirely mathtext, so the
    # math part carries \mathbf and the plain-text part is bolded by fontweight.
    # \log is written out as letters inside \mathbf so it bolds as upright bold text.
    plt.ylabel(r'$\mathbf{-\,log_{10}}$ (adjusted p-value)',
               fontsize=axis_title_text_size, labelpad=0.1, fontweight='bold')

    plt.title(title, fontsize=title_text_size, fontweight='bold')
    plt.xticks(fontsize=axis_tick_text_size)
    plt.yticks(fontsize=axis_tick_text_size)

    if not show_legend or legend_loc == "none":
        plt.legend([], frameon=False)
    else:
        plt.legend(title='Significance',
                   handletextpad=0.0, frameon=False,
                   borderpad=0, columnspacing=0.5, alignment="left",
                   labelspacing=0.0,
                   title_fontproperties={'weight': 'bold', 'size': legend_text_size},
                   fontsize=legend_text_size,
                   loc=legend_loc, bbox_to_anchor=legend_bbox_to_anchor)
    
    if xlims is not None:
        plt.xlim(xlims)
    if ylims is not None:
        plt.ylim(ylims)
    if savefig is not None:
        plt.savefig(savefig, format='pdf', dpi=300)
    # plt.show()

    # return lists only if requested ----------------------------------------
    if return_gene_lists:                                            # CHANGED
        down_reg = [g for g, _ in sorted(down_reg,
                                         key=lambda x: abs(x[1]), reverse=True)]
        up_reg   = [g for g, _ in sorted(up_reg,
                                         key=lambda x: abs(x[1]), reverse=True)]
        return down_reg, up_reg
    return None                                                      # CHANGED


def plot_volcano_with_colors(
    dge_results_df, 
    lfc_threshold=0.5, 
    pval_threshold=0.01,
    group=None, 
    pval_col='pvals_adj', 
    logfc_col='logfoldchanges', 
    name_col='names', 
    colorby='celltype',
    title='Volcano Plot',
    xlabel=None,
    all_text_size=None,
    text_size=6,
    title_text_size=6,
    legend_text_size=6,
    legend_loc=6,
    color_legend_bbox_to_anchor=(1.0, 0.7),
    shape_legend_bbox_to_anchor=(1.0, 0.3),
    axis_title_text_size=6,
    axis_tick_text_size=6,
    xlims=None, 
    ylims=None,
    savefig=None, 
    show_grid=False,
    point_alpha=0.5,
    point_size=10, 
    figsize=(10, 6),
    palette="tab10",
    rasterize=True,
    show_gene_labels=True,         # NEW: toggle gene‐name annotation
    label_lfc_threshold=None,      # NEW: logFC threshold for labels
    label_pval_threshold=None,     # NEW: p-value threshold for labels
    return_gene_lists=True,        # NEW: toggle return of gene lists
    show_color_legend=True,        # NEW: hide cell-type colour legend
    show_shape_legend=True,        # NEW: hide significance shape legend
    # -----------------------------------------------------------------------
):
    """
    Plot a volcano plot for differential expression analysis results with cell types as colors
    and significance as shapes.

    Returns
    -------
    dict or None
        Dictionary of significant genes per cell type if return_gene_lists is True, else None.
    """
    # ---------- NEW default label thresholds --------------------------------
    if label_lfc_threshold is None:
        label_lfc_threshold = lfc_threshold
    if label_pval_threshold is None:
        label_pval_threshold = pval_threshold
    # ------------------------------------------------------------------------

    if all_text_size is not None:  # CHANGED
        text_size = title_text_size = legend_text_size = \
                    axis_title_text_size = axis_tick_text_size = all_text_size
    
    df              = dge_results_df
    pvals           = df[pval_col]
    logfoldchanges  = df[logfc_col]
    gene_names      = df[name_col]
    celltypes       = df[colorby]

    sig_pval, sig_logfc = pval_threshold, lfc_threshold

    # color palette
    unique_celltypes = df[colorby].astype("category").cat.categories
    if isinstance(palette, str):
        cmap = plt.get_cmap(palette)
        colors = cmap.colors[:len(unique_celltypes)] \
                 if hasattr(cmap, 'colors') else [cmap(i) for i in range(len(unique_celltypes))]
        color_map = {ct: colors[i % len(colors)] for i, ct in enumerate(unique_celltypes)}
    elif isinstance(palette, dict):
        color_map = palette
    else:
        raise ValueError("Palette must be a colormap string or a dict mapping cell types to colors.")

    markers = {'significant_up': 'o',
               'significant_down': 's',
               'non_significant': 'x'}

    plt.figure(figsize=figsize)
    plt.grid(show_grid)
    plt.rcParams.update({'font.size': text_size})

    significant_genes_by_celltype = {}

    color_handles, shape_handles = [], []

    for celltype in unique_celltypes:
        subset = df[df[colorby] == celltype]
        pvals_s, lfc_s, genes_s = subset[pval_col], subset[logfc_col], subset[name_col]

        sig_up   = (pvals_s < sig_pval) & (lfc_s >  sig_logfc)
        sig_down = (pvals_s < sig_pval) & (lfc_s < -sig_logfc)
        nonsig   = ~(sig_up | sig_down)

        plt.scatter(lfc_s[sig_up],   -np.log10(pvals_s[sig_up]),
                    alpha=point_alpha, s=point_size,
                    marker=markers['significant_up'],   c=[color_map[celltype]],
                    label=f'{celltype} (Up)',   rasterized=rasterize)
        plt.scatter(lfc_s[sig_down], -np.log10(pvals_s[sig_down]),
                    alpha=point_alpha, s=point_size,
                    marker=markers['significant_down'], c=[color_map[celltype]],
                    label=f'{celltype} (Down)', rasterized=rasterize)
        plt.scatter(lfc_s[nonsig],   -np.log10(pvals_s[nonsig]),
                    alpha=point_alpha, s=point_size,
                    marker=markers['non_significant'],  c=[color_map[celltype]],
                    label=f'{celltype} (Non-significant)', rasterized=rasterize)

        if celltype not in [h.get_label() for h in color_handles]:
            color_handles.append(plt.Line2D([0], [0], lw=0, marker='o',
                                            color=color_map[celltype], label=celltype))

        # ---------- collect genes if needed ---------------------------------
        if show_gene_labels or return_gene_lists:  # NEW
            down, up = [], []
            mask = (pvals_s < label_pval_threshold) & \
                   (np.abs(lfc_s) > label_lfc_threshold)
            for g, pv, lfc in zip(genes_s[mask], pvals_s[mask], lfc_s[mask]):
                (down if lfc < -label_lfc_threshold else up).append((g, lfc))
            significant_genes_by_celltype[celltype] = {
                'less_than_thr'  : sorted(down, key=lambda x: abs(x[1]), reverse=True),
                'greater_than_thr': sorted(up,   key=lambda x: abs(x[1]), reverse=True)
            }
        # --------------------------------------------------------------------

    # legend handles for shapes
    for key, m in markers.items():
        shape_handles.append(plt.Line2D([0], [0], lw=0, marker=m,
                                        color='black', label=key.replace('_', ' ').capitalize()))

    plt.axhline(-np.log10(sig_pval), color='black', linestyle='--')
    plt.axvline( sig_logfc, color='black', linestyle='--')
    plt.axvline(-sig_logfc, color='black', linestyle='--')

    plt.tick_params(axis='both', which='major', pad=0.2, size=5)
    plt.xticks(fontsize=axis_tick_text_size)
    plt.yticks(fontsize=axis_tick_text_size)

    if xlims is not None:
        plt.xlim(xlims)
    if ylims is not None:
        plt.ylim(ylims)
        
    # ------------------- annotate genes -------------------------------------
    if show_gene_labels:  # NEW
        texts = []
        for ct, gd in significant_genes_by_celltype.items():
            for g, lfc in gd['less_than_thr'] + gd['greater_than_thr']:
                row = df[(df[name_col] == g) & (df[colorby] == ct)]
                if not row.empty:
                    texts.append(plt.text(lfc, -np.log10(row[pval_col].values[0]),
                                          g, fontsize=text_size, ha='center', va='center'))
        if texts:
            adjust_text(texts, arrowprops={'arrowstyle': '-', 'color': 'black'})
    # ------------------------------------------------------------------------

    plt.xlabel(xlabel if xlabel is not None else r'$log_2$ fold-change',
               fontsize=axis_title_text_size, labelpad=0.1)
    plt.ylabel(r'$-\,\log_{10}\,(\mathrm{adjusted}\;p\mathrm{-}value)$',
               fontsize=axis_title_text_size, labelpad=0.1)
    plt.title(title, fontsize=title_text_size, fontweight='bold')

    # Add separate legends (each can be hidden independently)
    color_leg = None
    if show_color_legend:
        color_leg = plt.legend(handles=color_handles, title="Cell Types", handletextpad=0.0, frameon=False, borderpad=0, columnspacing=0.2, alignment="left",
                       labelspacing=0.0, title_fontproperties={'weight': 'bold', 'size': legend_text_size}, fontsize=legend_text_size, loc = 6, bbox_to_anchor=color_legend_bbox_to_anchor)
    if show_shape_legend:
        if color_leg is not None:
            # Preserve the colour legend before adding the second one
            plt.gca().add_artist(color_leg)
        plt.legend(handles=shape_handles, title="Significance", handletextpad=0.0, frameon=False, borderpad=0, columnspacing=0.2, alignment="left",
                       labelspacing=0.0, title_fontproperties={'weight': 'bold', 'size': legend_text_size}, fontsize=legend_text_size, loc = 6, bbox_to_anchor=shape_legend_bbox_to_anchor)
    
    
    
    if savefig is not None:
        plt.savefig(savefig, format='pdf', dpi=300)
    # plt.show()

    # ------------------ return value control --------------------------------
    if return_gene_lists:          # NEW
        return significant_genes_by_celltype
    return None                    # NEW


def plot_volcano_with_celltypes_and_groups(
    dge_results_df,
    lfc_threshold=0.5,
    pval_threshold=0.01,
    pval_col='pvals_adj',
    logfc_col='logfoldchanges',
    name_col='names',
    celltype_col='celltype',
    order = None,
    group_col='group',
    title='Volcano Plot',
    xlabel=None,
    all_text_size=None,
    text_size=6,
    title_text_size=6,
    legend_text_size=6,
    color_legend_bbox_to_anchor=(1.0, 0.75),
    shape_legend_bbox_to_anchor=(1.0, 0.35),
    axis_title_text_size=6,
    axis_tick_text_size=6,
    xlims=None,
    ylims=None,
    savefig=None,
    show_grid=False,
    point_alpha=0.5,
    point_size=10,
    figsize=(10, 6),
    palette="tab10",
    rasterize=True,
    show_gene_labels=True,
    label_lfc_threshold=None,
    label_pval_threshold=None,
    return_gene_lists=True,
    gray_non_sig_points=False,
    # NEW: grey-out non-significant points
    show_color_legend=True,        # NEW: hide cell-type colour legend
    show_shape_legend=True,        # NEW: hide group shape legend
):
    """
    Volcano plot coloured by *cell type* and shaped by *group*.

    If ``gray_non_sig_points`` is True, any point that does **not** pass both
    the logFC and p-value thresholds is shown in grey, irrespective of cell type.
    """
    # ---------- threshold defaults -----------------------------------------
    if label_lfc_threshold is None:
        label_lfc_threshold = lfc_threshold
    if label_pval_threshold is None:
        label_pval_threshold = pval_threshold

    # ---------- unify text sizes -------------------------------------------
    if all_text_size is not None:
        text_size = title_text_size = legend_text_size = \
            axis_title_text_size = axis_tick_text_size = all_text_size

    df = dge_results_df.copy()
    sig_pval, sig_lfc = pval_threshold, lfc_threshold

    # ---------- colour map (cell types) ------------------------------------
    if order is None:
        celltypes = df[celltype_col].astype("category").cat.categories
    else:
        celltypes = order
        
    if isinstance(palette, str):
        cmap = plt.get_cmap(palette)
        base_colours = cmap.colors if hasattr(cmap, "colors") else [cmap(i) for i in range(len(celltypes))]
        colour_map = {ct: base_colours[i % len(base_colours)] for i, ct in enumerate(celltypes)}
    elif isinstance(palette, dict):
        colour_map = palette
    else:
        raise ValueError("palette must be a colormap name or a dict.")

    # ---------- shape map (groups) -----------------------------------------
    groups = df[group_col].unique()
    base_markers = ['o', 's', '^', 'v', 'D', 'P', '*', 'X', '<', '>']
    marker_map = {g: base_markers[i % len(base_markers)] for i, g in enumerate(groups)}

    # ---------- set up figure ----------------------------------------------
    plt.figure(figsize=figsize)
    plt.grid(show_grid)
    plt.rcParams.update({'font.size': text_size})

    colour_handles, shape_handles = [], []
    significant_genes_by_celltype = {}

    # ---------- plotting loop ----------------------------------------------
    for ct in celltypes:
        for gp in groups:
            sub = df[(df[celltype_col] == ct) & (df[group_col] == gp)]
            if sub.empty:
                continue

            sig_mask = (sub[pval_col] < sig_pval) & (abs(sub[logfc_col]) > sig_lfc)

            # --- points that meet BOTH thresholds ---------------------------
            plt.scatter(sub.loc[sig_mask, logfc_col],
                        -np.log10(sub.loc[sig_mask, pval_col]),
                        c=[colour_map[ct]],
                        marker=marker_map[gp],
                        alpha=point_alpha, s=point_size,
                        rasterized=rasterize)

            # --- optionally grey-out the rest -------------------------------  # NEW
            if gray_non_sig_points:
                plt.scatter(sub.loc[~sig_mask, logfc_col],
                            -np.log10(sub.loc[~sig_mask, pval_col]),
                            c=['darkgray'],
                            marker=marker_map[gp],
                            alpha=point_alpha, s=point_size,
                            rasterized=rasterize)
            else:  # keep original colour
                plt.scatter(sub.loc[~sig_mask, logfc_col],
                            -np.log10(sub.loc[~sig_mask, pval_col]),
                            c=[colour_map[ct]],
                            marker=marker_map[gp],
                            alpha=point_alpha, s=point_size,
                            rasterized=rasterize)
            # ----------------------------------------------------------------

            # legend handles (added only once per label)
            if ct not in [h.get_label() for h in colour_handles]:
                colour_handles.append(plt.Line2D([0], [0], lw=0, marker='o',
                                                 color=colour_map[ct], label=ct))
            if gp not in [h.get_label() for h in shape_handles]:
                shape_handles.append(plt.Line2D([0], [0], lw=0, marker=marker_map[gp],
                                                 color='black', label=gp))

            # gene list collection
            if show_gene_labels or return_gene_lists:
                down, up = [], []
                lab_mask = (sub[pval_col] < label_pval_threshold) & \
                           (abs(sub[logfc_col]) > label_lfc_threshold)
                for g, l in zip(sub.loc[lab_mask, name_col],
                                sub.loc[lab_mask, logfc_col]):
                    (down if l < -label_lfc_threshold else up).append((g, l))
                if ct not in significant_genes_by_celltype:
                    significant_genes_by_celltype[ct] = {'less_than_thr': [],
                                                         'greater_than_thr': []}
                significant_genes_by_celltype[ct]['less_than_thr'].extend(down)
                significant_genes_by_celltype[ct]['greater_than_thr'].extend(up)

    # ---------- reference lines --------------------------------------------
    plt.axhline(-np.log10(sig_pval), color='black', ls='--')
    plt.axvline( sig_lfc,  color='black', ls='--')
    plt.axvline(-sig_lfc,  color='black', ls='--')

    plt.xlabel(xlabel if xlabel is not None else r'$log_2$ fold-change',
               fontsize=axis_title_text_size, labelpad=0.1)
    plt.ylabel(r'$-\,\log_{10}\,(\mathrm{adjusted}\;p\mathrm{-}value)$', labelpad=0.1,
               fontsize=axis_title_text_size)
    plt.title(title, fontsize=title_text_size, weight='bold')
    plt.tick_params(axis='both', which='major', pad=0.2, size=5)
    plt.xticks(fontsize=axis_tick_text_size)
    plt.yticks(fontsize=axis_tick_text_size)

    if xlims is not None:
        plt.xlim(xlims)
    if ylims is not None:
        plt.ylim(ylims)
        
    # ---------- gene labels -------------------------------------------------
    if show_gene_labels:
        from adjustText import adjust_text
        texts = []
        for ct, lsts in significant_genes_by_celltype.items():
            for g, l in lsts['less_than_thr'] + lsts['greater_than_thr']:
                row = df[(df[name_col] == g) & (df[celltype_col] == ct)]
                if not row.empty:
                    texts.append(plt.text(l, -np.log10(row[pval_col].values[0]),
                                          g, fontsize=text_size, ha='center', va='center'))
        if texts:
            adjust_text(texts, arrowprops={'arrowstyle': '-', 'color': 'black'})

            
    # Add separate legends (each can be hidden independently)
    color_leg = None
    if show_color_legend:
        color_leg = plt.legend(handles=colour_handles, title="Cell Types", handletextpad=0.0, frameon=False, borderpad=0, columnspacing=0.5, alignment="left",
                       labelspacing=0.0, title_fontproperties={'weight': 'bold', 'size': legend_text_size}, fontsize=legend_text_size, loc = 6, bbox_to_anchor=color_legend_bbox_to_anchor)
    if show_shape_legend:
        if color_leg is not None:
            plt.gca().add_artist(color_leg)
        plt.legend(handles=shape_handles, title="Group", handletextpad=0.0, frameon=False, borderpad=0, columnspacing=0.5, alignment="left",
                       labelspacing=0.0, title_fontproperties={'weight': 'bold', 'size': legend_text_size}, fontsize=legend_text_size, loc = 6, bbox_to_anchor=shape_legend_bbox_to_anchor)

    if savefig is not None:
        plt.savefig(savefig, format='pdf', dpi=300)
    plt.show()

    # ---------- return gene lists ------------------------------------------
    if return_gene_lists:
        for ct in significant_genes_by_celltype:
            significant_genes_by_celltype[ct]['less_than_thr'] = [
                g for g, _ in sorted(significant_genes_by_celltype[ct]['less_than_thr'],
                                     key=lambda x: abs(x[1]), reverse=True)]
            significant_genes_by_celltype[ct]['greater_than_thr'] = [
                g for g, _ in sorted(significant_genes_by_celltype[ct]['greater_than_thr'],
                                     key=lambda x: abs(x[1]), reverse=True)]
        return significant_genes_by_celltype
    return None
