#!/bin/tcsh

#BSUB -J jobName[1-10]       #job name AND job array
#BSUB -n 1                   #number of cores: generally only need 1 to 5
#BSUB -R span[hosts=1]       #distribute across 1 node
#BSUB -W 36:00               #walltime limit: hours:minutes
#BSUB -o /share/.../trial_%I.out  # filepath for the output
#BSUB -e /share/.../error_%I.err  # filepath for error

module load R

Rscript r_file_name.r $LSB_JOBINDEX