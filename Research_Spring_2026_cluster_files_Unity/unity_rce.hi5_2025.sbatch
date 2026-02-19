#!/bin/bash

#SBATCH --job-name rce.hi5_2025
#SBATCH -o rce.hi5_2025_%a.out
#SBATCH -e rce.hi5_2025_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time 2-00:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nduan@umass.edu
#SBATCH --array=1-120

cd ~/Research/
module load r-rocker-ml-verse/4.5.1_cuda12.8.1+apptainer 
R --no-restore --no-save < rce.hi5_2025.R > rce.hi5_2025_"${SLURM_ARRAY_TASK_ID}".Rout 
