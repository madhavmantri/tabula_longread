#!/bin/bash
#SBATCH --job-name=pacbio_syncronize
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
# go to pacbio directory
cd /oak/stanford/groups/quake/mmantri/group.quake/tabula_longread/pacbio
# run snakemake
# snakemake --workflow-profile slurm --jobs 20 --rerun-triggers mtime 
snakemake -s syncronize.smk --profile slurm --jobs 100 --rerun-incomplete --rerun-triggers mtime --config recollapsed_directory=recollapsed
echo "Snakemake pipeline completed at $(date)"

