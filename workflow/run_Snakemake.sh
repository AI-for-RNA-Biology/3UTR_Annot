#!/bin/bash
#SBATCH --time 2-00:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12
#SBATCH --partition cpu
#SBATCH --account neurocirt

# Uncomment if you need to activate conda
#eval "$(/idiap/temp/${USER}/miniconda3/condabin//conda shell.bash hook)"
conda activate snakemake

# Uncomment if you need to do a SETSHELL
# . /idiap/resource/software/initfiles/shrc
# SETSHELL htk

snakemake --cores 1 --use-conda --rerun-trigger mtime --jobs 12
