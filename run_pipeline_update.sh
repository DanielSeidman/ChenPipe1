#!/bin/bash
#SBATCH -J sm
#SBATCH -o out
#SBATCH -e err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 2:00:00

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate snakemake

#this forces jobs with updated parameters to rerun
snakemake -R $(snakemake --list-params-changes) --snakefile workflow/Snakefile --profile ./profiles/slurm
