#!/bin/bash

#SBATCH --job-name rce_n10_p50
#SBATCH -o rce_n10_p50_%a.out
#SBATCH -e rce_n10_p50_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time 3-0:00
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nduan@umass.edu
#SBATCH --array=1-120

cd ~/Research/
module load R
R --no-restore --no-save < rce_n10_p50.R > rce_n10_p50_"${SLURM_ARRAY_TASK_ID}".Rout