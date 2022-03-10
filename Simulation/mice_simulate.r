library(gtools)

# w = as.numeric(Sys.getenv('LSB_JOBINDEX'))
w=2
print(w)
set.seed(w)

# Sample size
N <- 32

dt = 1/512

# Determining good lambda values for the states
# trial = mice_format[, c('state', 'delta', 'theta', 'alpha', 'beta')]
# is_lamba = trial[trial$state == 2, ]; colMeans(is_lamba)
# nrem_lamba = trial[trial$state == 3, ]; colMeans(nrem_lamba)
# rem_lamba = trial[trial$state == 4, ]; colMeans(rem_lamba)


pars = c(c(matrix(c(-2.26568339,  0.08766060, 
                    -1.22022878, -4.44888558, 
                    -1.56180104, -0.08262607, 
                    -2.20978996, 0.05404948,
                    -2.41222255, 0.10833734,
                    -2.26568339,  0.08766060, 
                    -1.22022878, -4.44888558, 
                    -1.56180104, -0.08262607, 
                    -2.20978996, 0.05404948,
                    -2.41222255, 0.10833734,
                    -2.20978996, 0.05404948,), ncol=2, byrow=T)),
         c(-3.444682, -3.850148, -4.543295,
           -3.218876, -1.321756, -3.624341,
           -3.624341, -1.321756, -3.218876,
           -4.543295, -3.850148, -3.444682),
         c( 0.5, 0.5),
         c(1, 1, 1, 1),
         c(1.25, 1, 0.5, 0),
         c(1.5, 0.8, 0.5, 0),
         c(1, 1, 0.5, 0))

pars = c(c(matrix(c(-2.26568339,  0.08766060,
                    -1.22022878, -4.44888558,
                    -1.56180104, -0.08262607,
                    -2.20978996,  0.05404948,
                    -2.41222255,  0.10833734,
                    -1.9       ,  0.07      ,
                    -1.0       ,  0.1       ,
                    -1.0       ,  0.1       ,
                    -1.0       ,  0.1       ,
                    -1.0       ,  0.1       ,
                    -1.0       ,  0.1       ), ncol=2, byrow=T)),
         c(-3.444682, -3.850148, -4.543295,
           -3.218876, -1.321756, -3.624341,
           -3.624341, -1.321756, -3.218876,
           -4.543295, -3.850148, -3.444682),
         c( -6.52842355, -6.15970066),
         c(1, 1, 1, 1),
         c(1.25, 1, 0.5, 0),
         c(1.5, 0.8, 0.5, 0),
         c(1, 1, 0.5, 0))

par_index = list( beta=1:22, misclass=23:34, pi_logit=35:36, 
                  l_delta = 37:40, l_theta=41:44, l_alpha=45:48, l_beta=49:52)

init_logit = c( 1, exp(pars[par_index$pi_logit][1]), 
                exp(pars[par_index$pi_logit][2]), 0)
initProbs = init_logit / sum(init_logit)

resp_fnc = matrix(c(1, exp(pars[par_index$misclass][1]), 
                    exp(pars[par_index$misclass][2]), exp(pars[par_index$misclass][3]),
                    exp(pars[par_index$misclass][4]), 1, 
                    exp(pars[par_index$misclass][5]), exp(pars[par_index$misclass][6]),
                    exp(pars[par_index$misclass][7]), exp(pars[par_index$misclass][8]), 
                    1, exp(pars[par_index$misclass][9]),
                    exp(pars[par_index$misclass][10]),exp(pars[par_index$misclass][11]),
                    exp(pars[par_index$misclass][12]), 1), ncol=4, byrow=TRUE)

errorMat = resp_fnc / rowSums(resp_fnc)

lambda_mat = matrix(c(pars[par_index$l_delta], pars[par_index$l_theta],
                      pars[par_index$l_alpha], pars[par_index$l_beta]),
                    nrow = 4, byrow = T)
colnames(lambda_mat) = c("delta", "theta", "alpha", "beta")
rownames(lambda_mat) = c("LIMBO", "IS", "NREM", "REM")

lambda_mat = exp(lambda_mat)

beta <- pars[par_index$beta]

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Collect information about the real data set.
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

load('../Data_format/mice_format_total.rda')

ptnum <- unique(mice_format$ptnum)
N_mice <- length(ptnum)
sleepTime <- NULL
NumObs <- rep(0,N_mice)
for(i in 1:N_mice){
  subject <- mice_format[mice_format$ptnum==ptnum[i],,drop=FALSE]
  
  # The number of observations for each subject.
  NumObs[i] <- nrow(subject)
  
  # The times between observations.
  sleepTime <- c( sleepTime, max(subject$t2))
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Fill states using trans. rate matrix and resp. fnc. estimated on the real data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

Q <- function(t,beta){
  
  betaMat = matrix(beta, ncol = 2, byrow = F) # determine the covariates
  
  q1  = exp( c(1,t) %*% betaMat[1,] )  # Transition from LIMBO to IS
  q2  = exp( c(1,t) %*% betaMat[2,] )  # Transition from LIMBO to NREM
  q3  = exp( c(1,t) %*% betaMat[3,] )  # Transition from IS    to LIMBO
  q4  = exp( c(1,t) %*% betaMat[4,] )  # Transition from IS    to NREM
  q5  = exp( c(1,t) %*% betaMat[5,] )  # Transition from IS    to REM
  q6  = exp( c(1,t) %*% betaMat[6,] )  # Transition from NREM  to LIMBO
  q7  = exp( c(1,t) %*% betaMat[7,] )  # Transition from NREM  to IS
  q8  = exp( c(1,t) %*% betaMat[8,] )  # Transition from NREM  to REM
  q9  = exp( c(1,t) %*% betaMat[9,] )  # Transition from REM   to LIMBO
  q10 = exp( c(1,t) %*% betaMat[10,] ) # Transition from REM   to IS
  q11 = exp( c(1,t) %*% betaMat[11,] ) # Transition from REM   to NREM
  
  qmat = matrix(c(  0,  q1,  q2,  0,
                    q3,   0,  q4, q5,
                    q6,  q7,   0, q8,
                    q9, q10, q11,  0),
                nrow = 4, byrow = T)
  diag(qmat) = -rowSums(qmat)
  
  return(qmat)
}

# Need to account for the proportions from the dirichlet function

rawData <- NULL
propDeaths_sim <- 0
NumObs_sim <- NULL
for(i in 1:N){
  
  print(i)
  
  # Sample for an initial state.
  trueState <- sample(1:4, size=1, prob=initProbs)
  
  # Sample how long the mouse is asleep
  timeAsleep = sample(sleepTime, size = 1)
  # Sample the remaining states until death.
  seconds = seq(0, timeAsleep, dt)
  timeCentScale = (seconds - mean(seconds)) / sd(seconds)
  
  s <- trueState
  resp <- matrix(nrow = length(seconds), ncol = 4)
  colnames(resp) = c('delta','theta','alpha','beta')
  
  for (ii in 1:(length(timeCentScale) - 1)){ # need to account for no absorbing states
    
    time1 = timeCentScale[ii]
    print(time1)
    # Infinitesimal transition rates.
    qmat <- Q(time1,beta)
    
    # Possible next states.
    moveToStates <- which(qmat[s,] > 0)
    
    # Sample the wait times before transition to each of the next possible states.
    waitTimes <- rexp( n=length(moveToStates), rate= qmat[s,moveToStates])
    
    # If any of the wait times are smaller than dt, then transition to the state with the minimum wait time.
    min_waitTime <- min(waitTimes)
    if(min_waitTime < dt){  s <- moveToStates[ which(waitTimes == min_waitTime) ]  }
    
    trueState <- c( trueState, s)
    
    temp_resp <- rdirichlet(1, lambda_mat[s, ])
    resp[ii, ] <-  temp_resp
    
  }
  
  # Sample inter-observation times from the cav data set.  Maximum of 20 seconds in study.
  visitTimes <- seconds
  
  # If first visit time occurs after death, then subject is NOT entered into the study.
  n_i <- length(visitTimes)
  state <- trueState

  ptnum <- rep(i,n_i)
  seconds <- visitTimes
  rawData <- rbind( rawData, data.frame(ptnum,timeCentScale,seconds,state, resp) )
  
  NumObs_sim <- c( NumObs_sim, n_i)
  
}

temp_df = data.frame("ptnum" = NA, "t1" = NA, "t2" = NA, "t1_cs" = NA, 
                     "t2_cs" = NA, "state" = NA, 'delta' = NA,'theta' = NA,
                     'alpha' = NA, 'beta' = NA)
rowNum = 1
# Split into 5 second intervals
for (j in unique(rawData$ptnum)) {
  print(j)
  
  subject_j = rawData[rawData$ptnum == j, ]
  
  index_seq = seq(0, nrow(subject_j), 2048)
  
  for(i in 1:(length(index_seq) - 1)) {
    print(i)
    
    s_ind = index_seq[i] + 1
    e_ind = index_seq[i+1]
    test = subject_j[s_ind:e_ind, ] # grabs data every 5 seconds
    
    t1 = head(test$seconds,1)
    t2 = tail(test$seconds,1)
    
    t1_cs = head(test$timeCentScale,1)
    t2_cs = tail(test$timeCentScale,1)

    s = ""
    
    for(jj in 1:length(unique(test$state))) {
      temp_s = unique(test$state)[jj]
      
      s = paste0(s, temp_s)
    }
    
    temp_df[rowNum,] = c(j, as.numeric(t1), as.numeric(t2), as.numeric(t1_cs), 
                         as.numeric(t2_cs), as.numeric(s), colMeans(test)[5:8])
    rowNum = rowNum + 1
  }
}

temp_df$state[temp_df$state > 4] = 99
print(length(which(temp_df$state == 99)) / nrow(temp_df)) # This is the proportion of 99's
print(length(which(mice_format$state == 99)) / nrow(mice_format)) # This is the proportion of 99's

rawData = temp_df

N <- length(unique(rawData$ptnum))

for(i in 1:nrow(rawData)){	
  if(rawData$state[i] != 99) {
    rawData$state[i] <- sample(1:4, size=1, prob=errorMat[rawData$state[i],])  
  }
}

obstrue <- rep(0,nrow(rawData))

rawData <- cbind(rawData,obstrue)

miceData = rawData

rownames(miceData) <- NULL

save(miceData, file=paste("Data_simulation/miceData", w, ".rda", sep=''))

# Transition frequencies for the simulated data set.
# nTrans_sim <- rep(0,6)
# for(i in 1:(nrow(miceData) - 1)){
#     if(miceData$state[i] == 1 & miceData$state[i+1] == 2) {nTrans_sim[1] = nTrans_sim[1] + 1}
#     else if (miceData$state[i] == 1 & miceData$state[i+1] == 3) {nTrans_sim[2] = nTrans_sim[2] + 1}
#     else if (miceData$state[i] == 2 & miceData$state[i+1] == 1) {nTrans_sim[3] = nTrans_sim[3] + 1}
#     else if (miceData$state[i] == 2 & miceData$state[i+1] == 3) {nTrans_sim[4] = nTrans_sim[4] + 1}
#     else if (miceData$state[i] == 3 & miceData$state[i+1] == 1) {nTrans_sim[5] = nTrans_sim[5] + 1}
#     else if (miceData$state[i] == 3 & miceData$state[i+1] == 2) {nTrans_sim[6] = nTrans_sim[6] + 1}
# }
# cat('Simulated data set sample size                         = ', N,'\n')
# cat('Simulated data set transition fequencies               = ', nTrans_sim / sum(nTrans_sim),'\n')
# cat('Simulated data set transition counts                   = ', nTrans_sim,'\n')
# cat('Simulated data set quantiles of number of observations = ','\n')
# print(quantile(NumObs_sim))
# print(sum(miceData$obstrue == 0))
