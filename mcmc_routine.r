library(mvtnorm, quietly=T)
library(foreach, quietly=T)
library(doParallel, quietly=T)

library(msm)
library(deSolve)
library(expm)
library(gtools)

Q <- function(t,x_ik,beta){

    betaMat = matrix(beta, ncol = 4, byrow = F) # determine the covariates
  
    q1  = exp( c(1,t,x_ik) %*% betaMat[1,] )  # Transition from AWAKE to IS
    q2  = exp( c(1,t,x_ik) %*% betaMat[2,] )  # Transition from AWAKE to NREM
    q3  = exp( c(1,t,x_ik) %*% betaMat[3,] )  # Transition from IS    to AWAKE
    q4  = exp( c(1,t,x_ik) %*% betaMat[4,] )  # Transition from IS    to NREM
    q5  = exp( c(1,t,x_ik) %*% betaMat[5,] )  # Transition from IS    to REM
    q6  = exp( c(1,t,x_ik) %*% betaMat[6,] )  # Transition from NREM  to AWAKE
    q7  = exp( c(1,t,x_ik) %*% betaMat[7,] )  # Transition from NREM  to IS
    q8  = exp( c(1,t,x_ik) %*% betaMat[8,] )  # Transition from NREM  to REM
    q9  = exp( c(1,t,x_ik) %*% betaMat[9,] )  # Transition from REM   to IS
    q10 = exp( c(1,t,x_ik) %*% betaMat[10,] ) # Transition from REM   to NREM
    
    qmat = matrix(c(  0,  q1,  q2,  0,
                     q3,   0,  q4, q5,
                     q6,  q7,   0, q8,
                      0,  q9, q10,  0),
                nrow = 4, byrow = T)
    diag(qmat) = -rowSums(qmat)

  return(qmat)
}

model_t <- function(t,p,parms) {

    qmat = Q(t, parms$x_ik, parms$b)

    pmat = matrix(c(  p[1],  p[2],  p[3],     0,
                      p[4],  p[5],  p[6],  p[7],
                      p[8],  p[9], p[10], p[11],
                         0, p[12], p[13], p[14],),
                nrow = 4, byrow = T)

    # Vectorizing the matrix multiplication row-wise
    dP = c(t(pmat %*% qmat))

    return(list(dP))

}

fn_log_post_continuous <- function(pars, prior_par, par_index, x, y_1, y_2, t, id) {

    # Order: Awake, IS, NREM, REM
    init_logit = c( 1, exp(pars[par_index$pi_logit][1]), 
                    exp(pars[par_index$pi_logit][2]), 0)
    init = init_logit / sum(init_logit)

    # the misclassification response function looks like the following
    # (0.94, 0.03, 0.02, 0.01)
    # (0.03, 0.75, 0.20, 0.02)
    # (0.02, 0.20, 0.75, 0.03)
    # (0.01, 0.02, 0.03, 0.94)
    resp_fnc = matrix(c(1, exp(pars[par_index$misclass][1]), 
            exp(pars[par_index$misclass][2]), exp(pars[par_index$misclass][3]),
            exp(pars[par_index$misclass][4]), 1, 
            exp(pars[par_index$misclass][5]), exp(pars[par_index$misclass][6]),
            exp(pars[par_index$misclass][7]), exp(pars[par_index$misclass][8]), 
            1, exp(pars[par_index$misclass][9]),
            exp(pars[par_index$misclass][10]),exp(pars[par_index$misclass][11]),
            exp(pars[par_index$misclass][12]), 1), ncol=4, byrow=TRUE)

    resp_fnc = resp_fnc / rowSums(resp_fnc)
    
    mu = pars[par_index$mu]
    sigma = pars[par_index$sigma]

    beta <- pars[par_index$beta]

    # Needs to be the same dimension as the dP in model_t() and represent the 
    #   identity matrix
    p_ic <- c(p1=1, p2=0, p3=0,
              p4=0, p5=1, p6=0, p7=0,
              p8=0, p9=0,p10=1,p11=0,
                   p12=0,p13=0,p14=1)

    log_total_val = foreach(i=unique(id), .combine='+', 
                            .export = c("model_t", "Q"), 
                            .packages = "deSolve") %dopar% {

        val = 1

        y_1_i = y_1[id == i]
        y_2_i = y_2[id == i]

        x_i = x[id == i,"sex",drop = F]

        t_i = t[id == i]

        d_1 = dnorm(y_2_i[1], mean = mu[1], sd = sigma)
        d_2 = dnorm(y_2_i[1], mean = mu[2], sd = sigma)
        d_3 = dnorm(y_2_i[1], mean = mu[3], sd = sigma)

        f_i = init %*% diag(c(d_1,d_2,d_3) * resp_fnc[, y_1_i[1]])
        log_norm = 0
        for(k in 2:length(t_i)) {

            out <- deSolve::ode(p_ic, times = t_i[(k-1):k], 
                                      func = model_t, 
                                      parms = list(b=beta, x_ik = x_i[k,]))

            P <- matrix(c(out[2,"p1"], out[2,"p2"], out[2,"p3"],
                        out[2,"p4"], out[2,"p5"], out[2,"p6"],
                        out[2,"p7"], out[2,"p8"], out[2,"p9"]), nrow = 3, byrow = T)

            d_1 = dnorm(y_2_i[k], mean = mu[1], sd = sigma)
            d_2 = dnorm(y_2_i[k], mean = mu[2], sd = sigma)
            d_3 = dnorm(y_2_i[k], mean = mu[3], sd = sigma)

            D_i = diag(c(d_1,d_2,d_3) * resp_fnc[, y_1_i[k]])

            val = f_i %*% P %*% D_i

            norm_val = sqrt(sum(val^2))
            f_i = val / norm_val
            log_norm = log_norm + log(norm_val)
        }

        return(log(sum(f_i)) + log_norm)
    }

    mean = prior_par$prior_mean
    sd = diag(prior_par$prior_sd)
    log_prior_dens = dmvnorm( x=pars, mean=mean, sigma=sd, log=T)

    return(log_total_val + log_prior_dens)

}



# -----------------------------------------------------------------------------
# The mcmc routine for samping the parameters
# -----------------------------------------------------------------------------
mcmc_routine = function( y_1, y_2, x, t, id, init_par, prior_par, par_index,
                         steps, burnin, n_cores){

  cl <- makeCluster(n_cores, outfile="")
  registerDoParallel(cl)

  pars = init_par
  n = length(y_1)
  n_par = length(pars)
  chain = matrix( 0, steps, n_par)

  group = list(c(par_index$beta, par_index$misclass, par_index$pi_logit,
                 par_index$mu, par_index$sigma))
  n_group = length(group)

  pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(group[[j]]))
  pscale = rep( .0001, n_group)
  accept = rep( 0, n_group)

  # Evaluate the log_post of the initial pars
  log_post_prev = fn_log_post_continuous( pars, prior_par, par_index, x, y_1, y_2, t, id)

  if(!is.finite(log_post_prev)){
    print("Infinite log-posterior; choose better initial parameters")
    break
  }

  # Begin the MCMC algorithm --------------------------------------------------
  chain[1,] = pars
  for(ttt in 2:steps){
    for(j in 1:n_group){

      # Propose an update
      ind_j = group[[j]]
      proposal = pars
      proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],sigma=pcov[[j]]*pscale[j])

      # Compute the log density for the proposal
      log_post = fn_log_post_continuous(proposal, prior_par, 
                                        par_index, x, y_1, y_2, 
                                        t, id)
      # print("Likelihood Evaluation:")
      # print(log_post)

      # Only propose valid parameters during the burnin period
      if(ttt < burnin){
        while(!is.finite(log_post)){
          print('bad proposal')
          proposal = pars
          proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],
                                     sigma=pcov[[j]]*pscale[j])
          log_post = fn_log_post_continuous(proposal, prior_par, 
                                            par_index, x, y_1, y_2, 
                                            t, id)
        }
      }

      # Evaluate the Metropolis-Hastings ratio
      if( log_post - log_post_prev > log(runif(1,0,1)) ){
        log_post_prev = log_post
        pars[ind_j] = proposal[ind_j]
        accept[j] = accept[j] +1
      }
      chain[ttt,ind_j] = pars[ind_j]
      print("Parameters Accepted")
      print(pars)

      # Proposal tuning scheme ------------------------------------------------
      if(ttt < burnin){
        # During the burnin period, update the proposal covariance in each step
        # to capture the relationships within the parameters vectors for each
        # transition.  This helps with mixing.
        if(ttt == 100)  pscale[j] = 1

        if(100 <= ttt & ttt <= 2000){
          temp_chain = chain[1:ttt,ind_j]
          pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])

        } else if(2000 < ttt){
          temp_chain = chain[(ttt-2000):ttt,ind_j]
          pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
        }
        if( sum( is.na(pcov[[j]]) ) > 0)  pcov[[j]] = diag( length(ind_j) )

        # Tune the proposal covariance for each transition to achieve
        # reasonable acceptance ratios.
        if(ttt %% 30 == 0){
          if(ttt %% 480 == 0){
            accept[j] = 0

          } else if( accept[j] / (ttt %% 480) < .4 ){ #change to 0.3
            pscale[j] = (.75^2)*pscale[j]

          } else if( accept[j] / (ttt %% 480) > .5 ){ #change to 0.4
            pscale[j] = (1.25^2)*pscale[j]
          }
        }
      }
      # -----------------------------------------------------------------------
    }
    # Restart the acceptance ratio at burnin.
    if(ttt == burnin)  accept = rep( 0, n_group)

    if(ttt%%1==0)  cat('--->',ttt,'\n')
  }
  # ---------------------------------------------------------------------------

  stopCluster(cl)
  print(accept/(steps-burnin))
  return(list( chain=chain[burnin:steps,], accept=accept/(steps-burnin),
               pscale=pscale))
}
# -----------------------------------------------------------------------------
