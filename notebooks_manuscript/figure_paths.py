#!/usr/bin/env python3
"""Central figure-output routing for the manuscript.

Every notebook saves via ``figpath()`` instead of a hard-coded ``figures/...``
string::

    plt.savefig(figpath("03_01_structural_category.pdf"), bbox_inches="tight", dpi=300)

``figpath`` looks the *basename* up in ``FIG_MAP``:

  - in the map  -> ``figures/<figure2|suppfigure7|...>/<name>``
  - not in map  -> ``figures/<name>``   (unchanged; exploratory figures stay put)

So reassigning a panel to a different manuscript figure is a one-line edit HERE,
not an edit to 68 notebooks. The target directory is created on demand.

FIG_MAP is regenerated from ``manuscript/figure_panel_map_REVIEW.csv`` by
``manuscript/apply_figure_map.py`` -- edit the CSV, not this file. Keys are BASENAMES.

.. warning::
   **A running kernel caches FIG_MAP.** ``figpath`` reads the module-level dict at
   call time, so a notebook that imported ``figure_paths`` *before* the map changed
   keeps writing to the old destination (usually ``figures/`` root). After
   ``apply_figure_map.py`` changes an assignment, either restart the kernel or::

       import importlib, figure_paths; importlib.reload(figure_paths)
       from figure_paths import figpath

   ``reload`` re-executes the module in place, so functions already imported pick up
   the new map. Re-running ``apply_figure_map.py`` also relocates anything that was
   written to the wrong place in the meantime.
"""
from pathlib import Path

# notebooks/ root, resolved from this file so it survives being imported from a
# subfolder (all notebooks chdir to notebooks/ via their .notebooks_root bootstrap,
# but be robust to callers that do not).
_ROOT = Path(__file__).resolve().parent
FIG_ROOT = _ROOT / "figures"

MAIN_FIGURES = 4
SUPP_FIGURES = 10

# basename -> subdirectory under figures/  ("figure1".."figure4", "suppfigure1".."suppfigure10")
# EMPTY = everything still lands in figures/ exactly as before. Populate after review.
FIG_MAP: dict[str, str] = {
    "02_01_pacbio_cell_counts_bubble_bysex.pdf": "figure1",
    "02_02_umap_broad_cell_class.pdf": "figure1",
    "03_01_structcat_proportion_vs_count.pdf": "figure2",
    "03_01_structural_category.pdf": "figure2",
    "03_02_AS_subcategory_catalog.pdf": "figure2",
    "04_01_novel_fraction_by_biotype.pdf": "figure2",
    "04_10_umap_hv_isoforms_celltype.pdf": "figure3",
    "04_11_case_study_FBXL13_expression.pdf": "figure3",
    "04_11_case_study_S100A16_expression.pdf": "figure3",
    "04_11_diu_de_venn.pdf": "figure3",
    "04_16_quadrant_subcategory_mix.pdf": "figure3",
    "05_01_variance_partition_diversity_renyi_boxplot.pdf": "figure3",
    "renyi_acorss_celltypes.pdf": "figure3",
    "renyi_entropy_mean_vs_variance.pdf": "figure3",
    "06_03_p16_vs_p14ARF_percell_scatter.pdf": "figure4",
    "06_08_fig4d_unconditional_depth_matched.pdf": "figure4",
    "06_10_selected_gene_percell_entropy_alltypes.pdf": "figure4",
    "06_10_volcano_depth_matched.pdf": "figure4",
    "cdkn2a_reference_isoforms.pdf": "figure4",
    "02_01_pacbio_umi_genes_molecules_per_cell_box_bysex.pdf": "suppfigure1",
    "04_05_scalar_vs_matrix_renyi.pdf": "suppfigure10",
    "04_09_renyi_vs_shannon_normshannon_tau.pdf": "suppfigure10",
    "04i_1minustau_vs_renyi_spearman_by_isoform.pdf": "suppfigure10",
    "05_01_resolution_ladder.pdf": "suppfigure10",
    "05_01_variance_partition_composition_dfadj.pdf": "suppfigure10",
    "05_01_variance_partition_composition_raw.pdf": "suppfigure10",
    "05_01_variance_partition_diversity_renyi_raw.pdf": "suppfigure10",
    "01_10_cdkn2a_pos_proportion_per_tissue_xgen_raw.pdf": "suppfigure11",
    "06_03_cdkn2a_molecule_enrichment_by_donor.pdf": "suppfigure11",
    "06_03_p16_p14ARF_venn.pdf": "suppfigure11",
    "06_07_fig5_codetection_forest.pdf": "suppfigure11",
    "06_08_supp11f_percelltype_entropy_matched.pdf": "suppfigure11",
    "06_09_supp11e_rarefaction_and_category.pdf": "suppfigure11",
    "06_11_CTSL_myeloid_leukocyte_percell_entropy_by_donor.pdf": "suppfigure12",
    "06_11_FN1_fibroblast_percell_entropy_by_donor.pdf": "suppfigure12",
    "06_11_top5_with_power_up.pdf": "suppfigure12",
    "06_12_pub_H_top5_diu_with_power.pdf": "suppfigure12",
    "02_01_sr_pacbio_barcode_overlap_bysex.pdf": "suppfigure2",
    "02_03_pacbio_T29_vs_T30_scatter.pdf": "suppfigure2",
    "02_03_shared_perdonor_scatter.pdf": "suppfigure2",
    "02_03_shared_persample_correlation.pdf": "suppfigure2",
    "03_01_coding_potential_by_structural.pdf": "suppfigure3",
    "03_01_exon_length_stacked_bars.pdf": "suppfigure3",
    "03_01_novel_across_tissues.pdf": "suppfigure3",
    "03_01_novel_isoform_sample_breadth_all_samples.pdf": "suppfigure3",
    "03_01_novel_isoform_sample_breadth_excluding_testis.pdf": "suppfigure3",
    "03_01_novel_vs_known_expression.pdf": "suppfigure3",
    "03_03_nmd_rate_by_structural_category.pdf": "suppfigure3",
    "04_01_novel_molecule_fraction_per_gene_distribution.pdf": "suppfigure3",
    "04_01_novel_molecules_per_gene_distribution.pdf": "suppfigure3",
    "04_01_ranked_isoforms_per_gene.pdf": "suppfigure3",
    "03_singlecell_functional_metrics_by_broad_cell_class.pdf": "suppfigure4",
    "03_singlecell_functional_metrics_by_tissue.pdf": "suppfigure4",
    "03_02_alt_splice_site_offsets.pdf": "suppfigure5",
    "03_02_novel_junction_sample_breadth_all_samples.pdf": "suppfigure5",
    "03_02_novel_junction_sample_breadth_excluding_testis.pdf": "suppfigure5",
    "03_02_novel_junctions_fractions_per_tissue.pdf": "suppfigure5",
    "03_02_splice_site_motifs.pdf": "suppfigure5",
    "03_04_TSS_TTS_offsets.pdf": "suppfigure5",
    "03_04_novel_TSS_CAGE_support.pdf": "suppfigure5",
    "03_04_novel_TTS_polyA_support.pdf": "suppfigure5",
    "03_04_per_gene_TSS_diversity.pdf": "suppfigure5",
    "03_06b_percell_LRpos_vs_LRneg_box.pdf": "suppfigure6",
    "03_06b_percell_LRpos_vs_LRneg_scatter.pdf": "suppfigure6",
    "03_06b_sr_junction_validation_rates_by_set.pdf": "suppfigure6",
    "03_06b_sr_support_depth_cells_allsamples.pdf": "suppfigure6",
    "03_06b_sr_support_depth_molecules_allsamples.pdf": "suppfigure6",
    "03_06b_validation_vs_molecule_count_reads_allsamples.pdf": "suppfigure6",
    "03_08_snaptron_support_level.pdf": "suppfigure6",
    "03_08_snaptron_validation_overall.pdf": "suppfigure6",
    "04_11_diu_breadth_isoform_vs_gene.pdf": "suppfigure7",
    "04_11_diu_de_quadrant_counts.pdf": "suppfigure7",
    "04_11_diu_vs_de_genes_per_cell_type.pdf": "suppfigure7",
    "04_11_diu_vs_de_hexbin.pdf": "suppfigure7",
    "04_17_fraction_vs_gene_rho_by_dominance_embedding.pdf": "suppfigure7",
    "04_17_fraction_vs_gene_rho_embedding.pdf": "suppfigure7",
    "04_17_knn_celltype_purity.pdf": "suppfigure7",
    "04_17_neighbour_overlap.pdf": "suppfigure7",
    "04_11_dm_diu_dotplot_delta_fraction_top_per_cell_class.pdf": "suppfigure8",
    "04_05_quadrant_example_genes_stackedbar_unnormalized.pdf": "suppfigure9",
    "04_16_quadrant_cds_cv.pdf": "suppfigure9",
    "04_16_quadrant_novel_fraction.pdf": "suppfigure9",
    "04_16_var_entropy_donor_reproducibility.pdf": "suppfigure9",
    "renyi_acorss_tissues.pdf": "suppfigure9",
    "renyi_entropy_distribution.pdf": "suppfigure9",
}


def figure_dirs() -> list[Path]:
    """The manuscript figure subdirectories (created if missing)."""
    names = ([f"figure{i}" for i in range(1, MAIN_FIGURES + 1)]
             + [f"suppfigure{i}" for i in range(1, SUPP_FIGURES + 1)])
    out = []
    for n in names:
        d = FIG_ROOT / n
        d.mkdir(parents=True, exist_ok=True)
        out.append(d)
    return out


def figpath(name: str) -> str:
    """Resolve a figure filename to its output path, creating the directory.

    ``name`` may be a bare basename or any path ending in one; only the basename
    is used for the FIG_MAP lookup, so existing ``figures/x.pdf`` strings can be
    passed through unchanged during migration.
    """
    base = Path(name).name
    sub = FIG_MAP.get(base)
    out = (FIG_ROOT / sub / base) if sub else (FIG_ROOT / base)
    out.parent.mkdir(parents=True, exist_ok=True)
    return str(out)


def unmapped(verbose: bool = True) -> list[str]:
    """Figure files sitting in figures/ root that FIG_MAP does not place."""
    if not FIG_ROOT.exists():
        return []
    loose = sorted(p.name for p in FIG_ROOT.iterdir()
                   if p.is_file() and p.name not in FIG_MAP)
    if verbose:
        print(f"{len(loose)} figure(s) in figures/ root not assigned to a manuscript figure")
    return loose


if __name__ == "__main__":
    dirs = figure_dirs()
    print(f"ensured {len(dirs)} figure subdirectories under {FIG_ROOT}")
    print(f"FIG_MAP entries: {len(FIG_MAP)}")
    unmapped()
