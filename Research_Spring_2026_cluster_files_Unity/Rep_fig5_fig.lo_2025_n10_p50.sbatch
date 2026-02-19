#!/bin/bash

#SBATCH --job-name Rep_fig5_fig.lo_2025_n10_p50
#SBATCH -o Rep_fig5_fig.lo_2025_n10_p50_%a.out
#SBATCH -e Rep_fig5_fig.lo_2025_n10_p50_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time 3-0:00
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nduan@umass.edu
#SBATCH --array=1-120

cd ~/Research/
module load R
R --no-restore --no-save < Rep_fig5_fig.lo_2025_n10_p50.R > Rep_fig5_fig.lo_2025_n10_p50_"${SLURM_ARRAY_TASK_ID}".Rout