#!/bin/bash
#SBATCH -J sm
#SBATCH -o out
#SBATCH -e err
#SBATCH -p RM-shared
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH --ntasks-per-node=128

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate snparcher
snakemake --snakefile workflow/Snakefile --profile ./profiles/slurm

