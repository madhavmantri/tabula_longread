#!/bin/bash
#SBATCH --job-name=pacbio_isoseq_xgen_sharded
#SBATCH --output=slurm_logs/process_%A.out
#SBATCH --error=slurm_logs/process_%A.err
#SBATCH --time=02-00:00:00
#SBATCH --mem=4G
#SBATCH --partition=quake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmantri@stanford.edu

# Sharded isoseq groupdedup for a single tissue. Plan → split → N×dedup → merge.
# WARNING: if the non-sharded groupdedup job for this tissue is still running,
# scancel it before launching this — both write to isoseq/dedup/{tissue}/.

# load environment
source /home/groups/quake/mmantri/miniconda3/etc/profile.d/conda.sh
conda activate pacbio

cd /oak/stanford/groups/quake/mmantri/group.quake/tabula_longread/pacbio_xgen

snakemake \
    -s groupdedup_sharded.smk \
    --profile slurm \
    --jobs 100 \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --config \
        hifi_locations_csv=/oak/stanford/groups/quake/mmantri/group.quake/tabula_longread/csvs/hifi_locations_xgen.csv \
        n_shards=100 \
        max_reads_per_barcode=10000 \
        percentile_cutoff=0

echo "Sharded groupdedup pipeline completed at $(date)"
