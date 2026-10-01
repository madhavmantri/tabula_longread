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

## Instructions for use

The pipelines are written for a SLURM cluster and were run on CentOS Linux 7.9.
Paths in the config files are absolute paths on the cluster where the work was
done, so repoint them before running.

1. **Create the environments.** The environment files pin the exact versions
   used for the manuscript.
   - `conda env create -f pacbio/pacbio.yml`: the pipeline (Snakemake, pbskera,
     lima, isoseq, pbmm2, pigeon, samtools). `pacbio_isoseq.sh` and
     `pacbio_syncronize.sh` activate it as `pacbio`.
   - `conda env create -f envs/sqanti3.yml`: the dependencies of SQANTI3 v6.0.1,
     which `pacbio_sqanti3.sh` activates as `sqanti3`. SQANTI3 itself is not a
     conda package: install v6.0.1 from its repository and point the
     `sqanti3_qc_script`, `sqanti3_filter_script` and `isoannotlite_script` keys
     in `pacbio/config.yaml` at it.
2. **List your samples.** Write a CSV with columns `tissue,hifi_dir`, one row per
   sample. `tissue` is the sample identifier used in every output path, and
   `hifi_dir` is the directory holding that sample's HiFi BAMs. Every
   `*.hifi_reads*.bam` in the directory is used, excluding `*.unassigned.bam`.
   If a directory holds BAMs for several samples, files ending `_{tissue}.bam`
   are matched to their own row. See `csvs/hifi_locations.csv` for the format.
   Write a second CSV, `sample_metadata.csv`, mapping each sample (`tube_id`)
   to its `donor` and `tissue`.
3. **Edit `pacbio/config.yaml`.**
   - `reference_genome`, `reference_gtf`: genome FASTA and GTF (the manuscript
     used the GENCODE GRCh38.p13 primary assembly and GENCODE v41).
   - `primer_fasta`: the Kinnex/MAS-Seq primer FASTA.
   - `hifi_locations_csv`, `sample_metadata_csv`: the two CSVs from step 2.
   - `work_directory`: where outputs are written.
   - `cage_peak`, `isoannotlite_script`, `sqanti3_qc_script`,
     `sqanti3_filter_script`: SQANTI3 data files and scripts.
   - Per-rule `mem`, `threads`, `time` and `partition`: set these for your
     cluster.
4. **Edit the launchers.** In `pacbio_isoseq.sh`, `pacbio_syncronize.sh` and
   `pacbio_sqanti3.sh`, update the `#SBATCH` lines, the conda activation line
   and the `cd` path.
5. **Run the three stages in order**, from `pacbio/`:
   ```bash
   sbatch pacbio_isoseq.sh       # per-sample processing: raw HiFi -> Seurat matrices
   sbatch pacbio_syncronize.sh   # cross-sample identifier harmonization (global re-collapse)
   sbatch pacbio_sqanti3.sh      # SQANTI3 QC, filtering and IsoAnnotLite
   ```
   Each launcher submits one Snakemake driver job, which submits the rule jobs.
   To see the planned jobs without running anything, run the same `snakemake`
   command with `-n`.
6. **Find the outputs.** Per sample, `pigeon/seurat/{tissue}/` holds the gene
   and isoform count matrices. After harmonization,
   `recollapsed/seurat/{tissue}/` holds the same matrices with atlas-wide
   transcript identifiers. These are the inputs to the notebooks in
   `notebooks_manuscript/01_atlas_construction/`.

For CDKN2A-targeted capture libraries, use `pacbio_xgen/` instead. Run its
three launchers in order: `pacbio_isoseq_xgen_prep.sh`,
`pacbio_isoseq_sharded.sh`, then `pacbio_isoseq_xgen_downstream.sh`.

## Reproducing the manuscript

Every notebook in `notebooks_manuscript/` is saved with its outputs, so each
reported number and figure panel can be checked without re-running anything.
To re-run notebooks, start from the deposited data:

- **Long-read processed data** on Figshare
  (https://doi.org/10.6084/m9.figshare.33411454): gene, pbid and ensemblid count
  objects, pbid-level pigeon and SQANTI3 classification, sample metadata, and
  per-cell CDKN2A capture calls.
- **Matched short-read data** from Tabula Sapiens, Gene Expression Omnibus
  accession GSE306755.

### Setup

```bash
git clone https://github.com/madhavmantri/tabula_longread.git
cd tabula_longread
conda env create -f envs/analysis.yml      # notebooks and build_from_figshare.py
conda activate analysis
mkdir -p pacbio/h5ads csvs
# put the three *_preprocessed.h5ad files from Figshare in pacbio/h5ads/
# and the other Figshare files (CSVs, GTF, BED) in csvs/
touch notebooks_manuscript/.notebooks_root
```

The lookup tables the notebooks read (`csvs/popv_to_ts_celltype_mapping.csv`,
`csvs/ts_celltypes.csv` and `csvs/transcripts_to_genes_with_biotypes.txt`)
ship with this repository.

`envs/popv.yml` is the environment for the cell-type annotation notebook
(`01_atlas_construction/01_07_popv_python3.11.ipynb`); every other notebook runs
in `envs/analysis.yml`.

Every notebook starts by walking up the directory tree to the
`.notebooks_root` marker and changing directory there. This makes relative
paths such as `./../pacbio/h5ads/` and `./../csvs/` resolve to the folders
created above. Git cannot track the marker, which is why the `touch` step is
needed.

### Building the derived objects

```bash
python notebooks_manuscript/build_from_figshare.py prepare
# run notebooks_manuscript/01_atlas_construction/01_09_compute_isoform_fraction.ipynb
python notebooks_manuscript/build_from_figshare.py cdkn2a
```

- **`prepare`** does two things:
  - It writes the annotated-only level. This is the ensemblid object restricted
    to features with an Ensembl transcript identifier, which gives the same
    145,051 features and 203,311 cells as in the manuscript.
  - It links the deposited classification tables to the file names the
    notebooks read.
  - It writes the per-cell popV predictions, which are stored on the objects.
  - It derives the ensemblid-level SQANTI3 table from the pbid-level one,
    keeping the first structure for each transcript. Features without an
    Ensembl match are identical to the manuscript table. For about 16% of
    annotated (ENST) features a different structure represents the
    transcript, because only the atlas structures are deposited, so
    `04_16` and `06_09` can differ slightly from the manuscript.
- **`01_09`** adds the isoform-fraction layer to all three isoform levels.
- **`cdkn2a`** writes the `*_xgen_cdkn2a.h5ad` objects used in
  `06_senescence/`. Each is its input object restricted to the 155,117 cells of
  the four capture donors, with the CDKN2A capture calls attached as `xgen_*`
  columns. The distinct-UMI counts are not deposited and are not rebuilt; only
  the whole-cell recovery comparison in `06_03` uses them.

The deposited SQANTI3 table renames `ORF_length` to `ORF_genomic_span`,
because that column holds the genomic span of the CDS, not its length.
Notebooks that still read `ORF_length` must use `CDS_length` instead.

### Short-read objects

Download the Tabula Sapiens objects from GEO GSE306755 and set
`DONOR_SR_PATHS` in `01_atlas_construction/01_06_merge_per_donor_shortread.ipynb`
to the downloaded files. `01_06` writes the merged and long-read-matched
short-read objects used by `02_01` and `02_03`.

### Running the notebooks

All notebooks can be run. The atlas-construction notebooks (`01_01`–`01_08`, `01_10`)
start from the pipeline outputs, so they need the raw sequencing data (see Data
availability in the manuscript). Every other notebook runs from the Figshare and GEO
data after the steps above, in folder order (`01` → `06`).
