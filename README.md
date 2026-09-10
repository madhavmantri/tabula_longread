# tabula_longread

Snakemake pipelines for PacBio MAS-Seq / Iso-Seq single-cell long-read
transcriptomics, as used to build a cross-tissue single-cell isoform atlas and a
CDKN2A-targeted capture dataset.

Two pipelines and the analysis notebooks behind the manuscript are published
here:

| Directory | Purpose |
|---|---|
| [`pacbio/`](pacbio) | Whole-transcriptome Iso-Seq: raw HiFi reads to per-sample transcript models, cross-sample identifier harmonization, SQANTI3 annotation and Seurat count matrices |
| [`pacbio_xgen/`](pacbio_xgen) | The same pipeline adapted for CDKN2A-targeted hybridization capture, where the barcode depth distribution is heavy-tailed enough that deduplication has to be sharded |
| [`notebooks_manuscript/`](notebooks_manuscript) | Downstream analysis: the notebooks that build the atlas objects from the pipeline outputs and produce every figure panel in the manuscript |

## Scope of this repository

**This repository contains code only** — the two pipelines and the analysis
notebooks. Sequencing data, intermediate BAMs, count matrices, h5ad objects,
rendered figures and manuscript sources are not included. Paths in the config
files and in the notebooks point at locations on the cluster where the work was
run and will need to be repointed before use.

## `pacbio/` — whole-transcriptome Iso-Seq

```
config.yaml            reference paths, barcode whitelist, per-rule SLURM resources
pacbio.yml             conda environment
isoseq.smk             main pipeline
syncronize.smk         cross-sample transcript identifier harmonization
sqanti3.smk            SQANTI3 QC, filtering and functional annotation
pacbio_isoseq.sh       launcher: isoseq.smk
pacbio_syncronize.sh   launcher: syncronize.smk
pacbio_sqanti3.sh      launcher: sqanti3.smk
```

`isoseq.smk` runs, per sample, the standard PacBio single-cell chain: `skera
split` to separate concatenated cDNA segments, `lima` to remove and orient
primers, `isoseq tag` / `refine` / `correct` for barcode and UMI handling,
`isoseq groupdedup` for deduplication, `pbmm2` alignment to the reference, and
`isoseq collapse` to build per-sample transcript models. Outputs are classified
with `pigeon` and written as Seurat-format matrices at gene and isoform level.

Because per-sample collapse assigns transcript identifiers independently, the
same structure receives a different ID in every sample. `syncronize.smk`
resolves this by pooling the per-sample collapsed sequences, re-collapsing them
globally and rewriting the original per-sample outputs with the shared
identifiers, so that a transcript means the same thing atlas-wide.

### `pacbio/scripts/`

```
sync_pacbio_ids.py             build and apply cross-sample identifier maps
recollapse_from_pooled_fasta.py  prefix per-sample FASTAs, then rewrite outputs
                                 after the global re-collapse
subsample_cell_barcode.py      cap reads per cell barcode before deduplication
saturation_rarefaction.py      rarefaction curves for sequencing saturation
prefix_and_pool_pacbio_ids.py  superseded by recollapse_from_pooled_fasta.py;
                               retained for reference
```

## `pacbio_xgen/` — CDKN2A-targeted capture

A self-contained copy of the pipeline with its own `config.yaml`, snakefiles and
scripts, so that capture-specific parameters cannot leak into the
whole-transcriptome run.

```
config.yaml                       capture-specific settings
isoseq.smk                        pipeline, adapted for the capture
groupdedup_sharded.smk            sharded deduplication
syncronize.smk                    identifier harmonization
pacbio_isoseq_xgen_prep.sh        phase A launcher: through subsampling and QC
pacbio_isoseq_sharded.sh          phase B launcher: sharded deduplication
pacbio_isoseq_xgen_downstream.sh  phase C launcher: alignment onward
pacbio_syncronize.sh              launcher: syncronize.smk
```

Targeted capture concentrates reads onto very few barcodes, and `isoseq
groupdedup` is effectively single-threaded, so deduplicating a capture library
in one pass does not finish in reasonable time. `groupdedup_sharded.smk` splits
the reads by cell barcode into bins of comparable size, deduplicates the bins in
parallel and merges the result. The pipeline is therefore run in three phases,
one launcher each, rather than as a single driver job.

### `pacbio_xgen/sharding_scripts/`

```
plan_dedup_shards.py       longest-processing-time bin-packing of barcodes
                           into shards of comparable cost
split_bam_by_cb_shards.py  scatter the BAM into per-shard files
```

The read cap used when planning shards must match the cap applied during
preparation, since the plan is built from pre-subsample barcode counts.

`pacbio_xgen/scripts/` mirrors `pacbio/scripts/`.

## Running

Both pipelines are Snakemake workflows submitted to SLURM. The launcher scripts
wrap the corresponding snakefile with the cluster settings; per-rule resources
are declared in `config.yaml` and resolved with a config-first lookup that falls
back to inline defaults.

```bash
conda env create -f pacbio/pacbio.yml
cd pacbio && ./pacbio_isoseq.sh
```

Reference genome, annotation, primer FASTA and barcode whitelist paths are all
set in `config.yaml`.

## `notebooks_manuscript/` — downstream analysis

The 39 notebooks on the manuscript critical path, together with the eight
helper modules that ship with them. A notebook is included if it produces a
figure panel used in the manuscript, or if it writes an object or table that
such a notebook reads; exploratory and QC notebooks from the working tree are
not published. Cell outputs are kept, so each notebook can be read as a record
of what was run without re-executing it.

| Directory | Notebooks | Contents |
|---|---|---|
| `01_atlas_construction/` | 10 | Per-sample Seurat matrices to merged, filtered, cell-type-labelled AnnData objects at four feature levels |
| `02_atlas_overview/` | 3 | Atlas composition, per-cell quality metrics, and comparison against matched short-read data |
| `03_isoform_splicing_landscape/` | 7 | Novel isoform catalogue, splicing events, coding and NMD status, UTR diversity, junction validation |
| `04_isoform_diversity_and_identity/` | 9 | Isoforms per gene, isoform diversity and specificity metrics, differential isoform usage across cell types |
| `05_celltype_vs_tissue/` | 1 | Variance partition of isoform usage between cell type and tissue |
| `06_senescence/` | 9 | CDKN2A isoform structure and the p16-positive senescence programme |

The helper modules sit in the tree root rather than beside the notebooks. Each
notebook opens by walking up the directory tree to a marker file and changing
directory to the root, which is what places the helpers on the import path and
makes the relative paths to the pipeline outputs resolve from any depth.

```
figure_paths.py                 maps a figure basename to its manuscript panel directory
add_classification.py           merge pigeon or SQANTI3 classification into an AnnData
isoform_fraction.py             per-gene isoform fraction layer
differential_isoform_fraction.py  differential isoform usage, including the
                                  Dirichlet-multinomial gene-level test
senescence_robustness.py        depth matching and robustness checks for the senescence analysis
pyVolcano.py                    volcano plots
TS_colorDict.py                 colour palettes for tissue, donor, compartment,
                                assay, sex and structural category
mm_process_adata_for_sags.py    leave-one-donor-out and per-donor differential
                                expression for senescence-associated genes
```

The notebooks read the count matrices and classification tables written by the
pipelines above, so they are not runnable from a clone alone; they are published
as the record of how the reported numbers and panels were produced.
