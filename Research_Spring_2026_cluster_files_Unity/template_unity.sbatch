#!/bin/bash
#SBATCH --job-name eels
#SBATCH -o Out/eels_%a.out
#SBATCH -e Out/eels_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-04:00:00
#SBATCH --mem-per-cpu=1g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maryclaregri@umass.edu
#SBATCH --array=1-1000

cd /work/pi_maryclaregri_umass_edu/Stopping/Code

module load r-rocker-ml-verse/4.5.1_cuda12.8.1+apptainer 

R --no-restore --no-save < eels_mcmc.R > Out/eels_mcmc_"${SLURM_ARRAY_TASK_ID}".Rout