#!/bin/bash

#SBATCH -o Model_out/trial%a.out
#SBATCH --array=1-50
#SBATCH --nodes=1
#SBATCH --ntasks=17
#SBATCH --account=jantonelli
#SBATCH --qos=jantonelli-b
#SBATCH --job-name=contResp
#SBATCH --time=01:00:00
#SBATCH -t 4000
#SBATCH --mem=5gb

export OMP_NUM_THREADS=17

module load R/4.0

R CMD BATCH --no-save mcmc_runfile_new.r Model_out/trial${SLURM_ARRAY_TASK_ID}.Rout
