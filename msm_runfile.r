library("msm")

# seedInd = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
seedInd = 5
set.seed(seedInd)

load('Data_format/mice_format_sub_total_split.rda')
mice_format = mice_format[!(mice_format$ptnum %in% 34:59), ]

sub_ind = c(101,  77,   3,  97,  70,  84, 113,  82,  30,  88,   2,  24,   1,  61,  27,  74,  99, 104,  32,
             73, 110,  15,  22,  90,  98,  75,  80,  11,  16,  14,  21,  78, 102,  91, 103,  96,  13,  20,
             12,  29,   5,  67,  68,  85,   8,  10,  60, 106,  92,  87)
 
mice_format = mice_format[(mice_format$ptnum %in% sub_ind), ]

obstrue <- rep(0,nrow(mice_format))
obstrue[which(mice_format$state == 99)] = 1
mice_format = cbind(mice_format, obstrue)

qmat <- matrix(c(0,     exp(-2),exp(-2),
                exp(-2),      0,exp(-2),
                exp(-2),exp(-2),      0), ncol=3, byrow=TRUE)
dimnames(qmat) <- list( c('IS','NREM','REM'), c('IS','NREM','REM'))

#----------------------------------------------------------------------------------------------------------------
# Run the msm implementation ------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

# emat = matrix(c(    1,   exp(-3), exp(-3),
#                 exp(-3),       1, exp(-3),
#                 exp(-3), exp(-3),       1), ncol=3, byrow=TRUE)
# emat = emat / rowSums(emat)
emat = diag(3)
dimnames(emat) <- list( c('IS','NREM','REM'), c('IS','NREM','REM'))

Output_msm <- msm(state ~ t1, subject=ptnum, data=mice_format, qmatrix=qmat, covariates= ~ 1 + t1, 
                center=FALSE, covinits=list(t1=rep(0,11)), obstrue=obstrue, 
                ematrix=emat, initprobs=c(1, exp(-0.523), exp(-0.7265), 0), est.initprobs=TRUE, deathexact=FALSE, 
                censor=99, censor.states=1:4, method='BFGS', control=list(fnscale=4000, maxit=10000))   

save(Output_msm,file=paste0('Model_out/msm/Output_msm_total_new',seedInd,'_15sec.rda'))
