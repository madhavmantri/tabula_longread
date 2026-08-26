#!/bin/bash
#SBATCH --job-name=pacbio_isoseq_xgen_downstream
#SBATCH --output=slurm_logs/process_%A.out
#SBATCH --error=slurm_logs/process_%A.err
#SBATCH --time=07-00:00:00
#SBATCH --mem=4G
#SBATCH --partition=quake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmantri@stanford.edu

# Phase B of the xgen pipeline: run every rule after isoseq_groupdedup
# (pbmm2_align → isoseq_collapse → pigeon_* → seurat). Assumes the sharded
# pipeline (pacbio_isoseq_sharded.sh) has already produced the merged
# isoseq/dedup/{tissue}/{tissue}.dedup.{bam,fasta} at the standard pipeline
# path. Snakemake skips groupdedup when those outputs exist and are newer
# than subsample.bam (the default `--rerun-triggers mtime` rule).
#
# WARNING: if you re-run pacbio_isoseq_xgen_prep.sh after sharded dedup
# finishes, the new subsample.bam mtime will be > dedup.bam mtime and
# snakemake will try to re-trigger isoseq_groupdedup. Re-run the sharded
# pipeline first in that case.

source /home/groups/quake/mmantri/miniconda3/etc/profile.d/conda.sh
conda activate pacbio

cd /oak/stanford/groups/quake/mmantri/group.quake/tabula_longread/pacbio_xgen

snakemake \
    -s isoseq.smk \
    --profile slurm \
    --jobs 100 \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --config \
hifi_locations_csv=/oak/stanford/groups/quake/mmantri/group.quake/tabula_longread/csvs/hifi_locations_xgen.csv \
        max_reads_per_barcode=10000 \
        percentile_cutoff=0

echo "xgen downstream pipeline completed at $(date)"
