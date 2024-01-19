#!/bin/bash
#SBATCH -J sm
#SBATCH -o out
#SBATCH -e err
#SBATCH -p RM
#SBATCH -n 1
#SBATCH -t 00:01:00

echo $(conda info --base)
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate snparcher
conda list
#snakemake --snakefile workflow/Snakefile --profile ./profiles/slurm
