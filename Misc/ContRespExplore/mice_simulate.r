library(msm)

num_iter = as.numeric(Sys.getenv('LSB_JOBINDEX'))
set.seed(num_iter)

# Set the sample size.  Note that the cav data set has 622 subjects.
N <- 2000
# Choose the discretization of time.
dt <- 1/365


# I need to change these paramters and parameter indices to account for the
# greater number of transition intensity estimates

trueValues <- c(c(matrix(c(-2.26568339,  0.08766060, -0.49991746,
                           -1.22022878, -4.44888558, -0.82779213,
                           -1.56180104, -0.08262607,  0.73838829,
                           -2.20978996,  0.05404948, -1.83682627,
                           -2.41222255,  0.10833734,  1.63135439,
                           -1.9       ,  0.07      , -1.1  ), ncol=3, byrow=T)),
                c(  -5.73343061, -2.140066, -2.833213, -2.12144526),
                c( -6.52842355, -6.15970066),
                c(10, 20, 30),
                1)

par_index = list( beta=1:18, misclass=19:22, pi_logit=23:24, mu = 25:27, sigma = 28)

betaMat <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)


errorMat_temp = matrix(c(1, exp(trueValues[par_index$misclass][1]), 0,
                         exp(trueValues[par_index$misclass][2]), 1, exp(trueValues[par_index$misclass][3]),
                         0, exp(trueValues[par_index$misclass][4]), 1), ncol=3, byrow=TRUE)

errorMat = errorMat_temp / rowSums(errorMat_temp)


initProbs_temp = c( 1, exp(trueValues[par_index$pi_logit][1]), exp(trueValues[par_index$pi_logit][2]), 0)
initProbs = initProbs_temp / sum(initProbs_temp)


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Collect information about the real data set.
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

ptnum <- unique(cav$PTNUM)
N_cav <- length(ptnum)
interObsTime <- NULL
propMale <- 0
propDeaths <- 0
NumObs <- rep(0,N_cav)
for(i in 1:N_cav){
  subject <- cav[cav$PTNUM==ptnum[i],,drop=FALSE]

  # The number of observations for each subject.
  NumObs[i] <- nrow(subject)

  # The times between observations.
  if(!(4 %in% subject$state)){  interObsTime <- c( interObsTime, round( diff(subject$years), 3))  }

  # Determine whether the subject is male.
  propMale <- propMale + as.integer(subject$sex[1]==1)

  # Determine whether the subject's death was observed.
  propDeaths <- propDeaths + as.integer(4 %in% subject$state)
}
propMale <- propMale / N_cav
propDeaths <- propDeaths / N_cav

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Fill in the states using the transition rate matrix and error matrix estimated on the real data set.
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

Q <- function(time,sex,betaMat){

  q1  = exp( c(1,time,sex) %*% betaMat[1,] )  # Transition from IS to NREM
  q2  = exp( c(1,time,sex) %*% betaMat[2,] )  # Transition from IS to REM
  q3  = exp( c(1,time,sex) %*% betaMat[3,] )  # Transition from NREM to IS
  q4  = exp( c(1,time,sex) %*% betaMat[4,] )  # Transition from NREM to REM
  q5  = exp( c(1,time,sex) %*% betaMat[5,] )  # Transition from REM to IS
  q6  = exp( c(1,time,sex) %*% betaMat[6,] )  # Transition from REM to NREM

  qmat = matrix(c( 0, q1,q2,
                  q3,  0,q4,
                  q5, q6, 0),nrow=3,byrow=TRUE)
  diag(qmat) = -rowSums(qmat)

  return(qmat)
}


rawData <- NULL
propDeaths_sim <- 0
NumObs_sim <- NULL
for(i in 1:N){

  print(i)

  # Sample the gender, as proportional to the cav data set.
  sex <- as.integer(runif(1,0,1) < propMale)

  # Sample for an initial state.
  trueState <- sample(1:4, size=1, prob=initProbs)

  # Sample the remaining states until death.
  years <- 0
  time1 <- 0
  s <- trueState

  while(time1 < 15){ # need to account for no absorbing states

    # Infinitesimal transition rates.
    qmat <- Q(time1,sex,betaMat)

    # Possible next states.
    moveToStates <- which(qmat[s,] > 0)

    # Sample the wait times before transition to each of the next possible states.
    waitTimes <- rexp( n=length(moveToStates), rate= qmat[s,moveToStates])

    # If any of the wait times are smaller than dt, then transition to the state with the minimum wait time.
    min_waitTime <- min(waitTimes)
    if(min_waitTime < dt){  s <- moveToStates[ which(waitTimes == min_waitTime) ]  }

    time1 <- time1 + dt

    years <- c( years, time1)
    trueState <- c( trueState, s)

  }
  timeOfDeath <- tail(years,1)

  # Sample inter-observation times from the cav data set.  Maximum of 20 years in study.
  visitTimes <- NULL
  time2 <- 0

  while(time2 < min( 20, timeOfDeath)){

    visitTimes <- c( visitTimes, time2)
    time2 <- time2 + sample( interObsTime, size=1)
  }

  # If first visit time occurs after death, then subject is NOT entered into the study.
  if( !is.null(visitTimes) ){

    # If death occured before the study ended, then record the time of death.
    if( timeOfDeath < 20 ){  visitTimes <- c( visitTimes, timeOfDeath) }


    n_i <- length(visitTimes)
    state <- NULL
    resp <- NULL
    for(k in 1:n_i){
        temp_state = tail( trueState[ years <= visitTimes[k] ], 1)
        state <- c( state, temp_state)
        temp_resp <- rnorm(1, trueValues[par_index$mu][temp_state],
                            trueValues[par_index$sigma])
        resp <-  c( resp, temp_resp)
    }

    ptnum <- rep(i,n_i)
    years <- visitTimes
    rawData <- rbind( rawData, data.frame(ptnum,years,sex,state, resp) )

    NumObs_sim <- c( NumObs_sim, n_i)
  }

}

colnames(rawData) <- c('ptnum','years','sex','state','cont_resp')
N <- length(unique(rawData$ptnum))
propDeaths_sim <- propDeaths_sim / N


# Add noise to the states.
for(i in 1:nrow(rawData)){	rawData$state[i] <- sample(1:3, size=1, prob=errorMat[rawData$state[i],])  }

obstrue <- rep(0,nrow(rawData))

hold <- cbind(rawData,obstrue)
hold <- hold[,c('ptnum','years','sex','state','cont_resp','obstrue')]

tempRow <- rep(0,ncol(hold))
names(tempRow) <- c('ptnum','years','sex','state','cont_resp','obstrue')

miceData = hold

colnames(miceData) <- c('ptnum','years','sex','state','cont_resp','obstrue')
rownames(miceData) <- NULL

save(miceData, file=paste("DataOut/Continuous/miceData", num_iter, ".rda", sep=''))

# Transition frequencies for the simulated data set.
nTrans_sim <- rep(0,6)
for(i in 1:(nrow(miceData) - 1)){
    if(miceData$state[i] == 1 & miceData$state[i+1] == 2) {nTrans_sim[1] = nTrans_sim[1] + 1}
    else if (miceData$state[i] == 1 & miceData$state[i+1] == 3) {nTrans_sim[2] = nTrans_sim[2] + 1}
    else if (miceData$state[i] == 2 & miceData$state[i+1] == 1) {nTrans_sim[3] = nTrans_sim[3] + 1}
    else if (miceData$state[i] == 2 & miceData$state[i+1] == 3) {nTrans_sim[4] = nTrans_sim[4] + 1}
    else if (miceData$state[i] == 3 & miceData$state[i+1] == 1) {nTrans_sim[5] = nTrans_sim[5] + 1}
    else if (miceData$state[i] == 3 & miceData$state[i+1] == 2) {nTrans_sim[6] = nTrans_sim[6] + 1}
}
cat('Simulated data set sample size                         = ', N,'\n')
cat('Simulated data set transition fequencies               = ', nTrans_sim / sum(nTrans_sim),'\n')
cat('Simulated data set transition counts                   = ', nTrans_sim,'\n')
cat('Simulated data set quantiles of number of observations = ','\n')
print(quantile(NumObs_sim))
print(sum(miceData$obstrue == 0))
