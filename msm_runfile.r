library("msm")

seedInd = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(seedInd)

load("Data_format/mice_format_total.rda")

obstrue <- rep(0,nrow(mice_format))
obstrue[which(mice_format$state == 99)] = 1
mice_format = cbind(mice_format, obstrue)

qmat <- matrix(c(0,     exp(-2),exp(-2),      0,
                exp(-2),      0,exp(-2),exp(-2),
                exp(-2),exp(-2),      0,exp(-2),
                exp(-2),exp(-2),exp(-2),      0), ncol=4, byrow=TRUE)
dimnames(qmat) <- list( c('LIMBO', 'IS','NREM','REM'), c('LIMBO', 'IS','NREM','REM'))

#----------------------------------------------------------------------------------------------------------------
# Run the msm implementation ------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

emat = matrix(c(    1,   exp(-3), exp(-3),   exp(-3),
                exp(-3),       1, exp(-3),   exp(-3),
                exp(-3), exp(-3),       1,   exp(-3),
                exp(-3), exp(-3), exp(-3),         1), ncol=4, byrow=TRUE)
emat = emat / rowSums(emat)
dimnames(emat) <- list( c('LIMBO', 'IS','NREM','REM'), c('LIMBO', 'IS','NREM','REM'))

Output_msm <- msm(state ~ t1, subject=ptnum, data=mice_format, qmatrix=qmat, covariates= ~ 1 + t1, 
                center=FALSE, covinits=list(t1=c(0,0,0,0,0)), obstrue=obstrue, 
                ematrix=emat, initprobs=c(1, exp(-0.523), exp(-0.7265), 0), est.initprobs=TRUE, deathexact=FALSE, 
                censor=99, censor.states=1:4, method='BFGS', control=list(fnscale=4000, maxit=10000))   

save(Output_msm,file=paste0('Model_out/msm/Output_msm',seedInd,'.rda'))
