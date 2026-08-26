#!/bin/bash
#SBATCH --job-name=pacbio_xgen_syncronize
#SBATCH --output=slurm_logs/process_%A.out
#SBATCH --error=slurm_logs/process_%A.err
#SBATCH --time=07-00:00:00
#SBATCH --mem=4G
#SBATCH --partition=quake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mmantri@stanford.edu

# load environment
source /home/groups/quake/mmantri/miniconda3/etc/profile.d/conda.sh
conda activate pacbio
# go to pacbio_xgen directory (self-contained copy: own config.yaml, scripts/, isoseq.smk)
cd /oak/stanford/groups/quake/mmantri/group.quake/tabula_longread/pacbio_xgen
# run snakemake — override hifi_locations_csv to the xGen tissue list (config.yaml
# defaults to the main hifi_locations.csv, same as the other xgen launchers)
snakemake -s syncronize.smk --profile slurm --jobs 100 --rerun-incomplete --rerun-triggers mtime --config \
    hifi_locations_csv=/oak/stanford/groups/quake/mmantri/group.quake/tabula_longread/csvs/hifi_locations_xgen.csv \
    recollapsed_directory=recollapsed
echo "Snakemake pipeline completed at $(date)"
