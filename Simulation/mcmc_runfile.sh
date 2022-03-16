#!/bin/tcsh

#BSUB -J miceModel[1-50]  #job name AND job array
#BSUB -n 17                  #number of cores
#BSUB -R span[hosts=1]       #distribute across 1 node
#BSUB -W 36:00               #walltime limit: hh:mm
#BSUB -o /share/hmmrs/ebkendal/HMM_Solver_Mice/Simulation/Model_out/a_out/trial_%I.out
#BSUB -e /share/hmmrs/ebkendal/HMM_Solver_Mice/Simulation/Model_out/a_out/error_%I.err  #error - %J is the job-id %I is the job-array index

module load R

Rscript mcmc_runfile.r $LSB_JOBINDEX