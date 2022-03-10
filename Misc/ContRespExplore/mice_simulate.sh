#!/bin/tcsh

#BSUB -J contResp[1-50]  #job name AND job array
#BSUB -n 5                   #number of cores
#BSUB -R span[hosts=1]       #distribute across 1 node
#BSUB -W 36:00               #walltime limit: hh:mm
#BSUB -o /share/hmmrs/ebkendal/HMM_Solver/Mice/DataOut/a_out/trial_%I.out
#BSUB -e /share/hmmrs/ebkendal/HMM_Solver/Mice/DataOut/a_out/error_%I.err  #error - %J is the job-id %I is the job-array index

module load R

Rscript mice_simulate.r $LSB_JOBINDEX
