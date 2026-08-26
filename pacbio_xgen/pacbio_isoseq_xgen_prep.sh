#!/bin/bash
#SBATCH --job-name=pacbio_isoseq_xgen_prep
#SBATCH --output=slurm_logs/process_%A.out
#SBATCH --error=slurm_logs/process_%A.err
#SBATCH --time=07-00:00:00
#SBATCH --mem=4G
#SBATCH --partition=quake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmantri@stanford.edu

# Phase A of the xgen pipeline: run every rule up to (and including)
# subsample_cell_barcode, isoseq_bcstats, and saturation_rarefaction.
# Stops before isoseq_groupdedup — that step is handled by
# pacbio_isoseq_sharded.sh. After sharded dedup finishes, run
# pacbio_isoseq_xgen_downstream.sh to complete the pipeline.

source /home/groups/quake/mmantri/miniconda3/etc/profile.d/conda.sh
conda activate pacbio

cd /oak/stanford/groups/quake/mmantri/group.quake/tabula_longread/pacbio_xgen

# --forcerun isoseq_correct isoseq_bcstats: percentile_cutoff is overridden below,
# but --rerun-triggers mtime won't detect a param-only change. Forcing these two
# rules makes the new percentile actually take effect. The cascade via mtime
# re-runs sort → bcstats/subsample/saturation, all bounded by --until.
# --forcerun isoseq_correct isoseq_bcstats \
# --until subsample_cell_barcode isoseq_bcstats saturation_rarefaction \

snakemake \
    -s isoseq.smk \
    --profile slurm \
    --jobs 100 \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --until subsample_cell_barcode isoseq_bcstats saturation_rarefaction \
    --forcerun subsample_cell_barcode \
    --config hifi_locations_csv=/oak/stanford/groups/quake/mmantri/group.quake/tabula_longread/csvs/hifi_locations_xgen.csv \
        percentile_cutoff=0 \
        max_reads_per_barcode=10000
echo "xgen prep pipeline completed at $(date)"
