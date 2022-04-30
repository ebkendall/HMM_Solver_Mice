#!/bin/tcsh

#BSUB -J mice[1-15]       #job name AND job array
#BSUB -n 16                   #number of cores: generally only need 1 to 5
#BSUB -R span[hosts=1]       #distribute across 1 node
#BSUB -W 48:00               #walltime limit: hours:minutes
#BSUB -o /share/hmmrs/ebkendal/HMM_Solver_Mice/FINAL/Model_out/a_out/trial1_%I.out  # filepath for the output
#BSUB -e /share/hmmrs/ebkendal/HMM_Solver_Mice/FINAL/Model_out/a_out/error1_%I.err  # filepath for error

module load R

Rscript mcmc_runfile.r $LSB_JOBINDEX
