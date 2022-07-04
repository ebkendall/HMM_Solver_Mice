library(mvtnorm, quietly=T)
library(foreach, quietly=T)
library(doParallel, quietly=T)

library(msm)
library(deSolve)
library(expm)
library(gtools)

Q <- function(t,beta){

    betaMat = matrix(beta, ncol = 2, byrow = F) # determine the covariates
  
    q1  = exp( c(1,t) %*% betaMat[1,] )  # Transition from IS    to NREM
    q2  = exp( c(1,t) %*% betaMat[2,] )  # Transition from IS    to REM
    q3  = exp( c(1,t) %*% betaMat[3,] )  # Transition from IS    to LIMBO
    q4  = exp( c(1,t) %*% betaMat[4,] )  # Transition from NREM  to IS
    q5  = exp( c(1,t) %*% betaMat[5,] )  # Transition from NREM  to REM
    q6  = exp( c(1,t) %*% betaMat[6,] )  # Transition from NREM  to LIMBO
    q7  = exp( c(1,t) %*% betaMat[7,] )  # Transition from REM   to IS
    q8  = exp( c(1,t) %*% betaMat[8,] )  # Transition from REM   to NREM
    q9  = exp( c(1,t) %*% betaMat[9,] )  # Transition from REM   to LIMBO
    q10 = exp( c(1,t) %*% betaMat[10,] ) # Transition from LIMBO to IS
    q11 = exp( c(1,t) %*% betaMat[11,] ) # Transition from LIMBO to NREM
    q12 = exp( c(1,t) %*% betaMat[12,] ) # Transition from LIMBO to REM
    
    qmat = matrix(c(  0,  q1,  q2, q3,
                     q4,   0,  q5, q6,
                     q7,  q8,   0, q9,
                    q10, q11, q12,  0),
                nrow = 4, byrow = T)
    diag(qmat) = -rowSums(qmat)

  return(qmat)
}

fn_log_post_continuous <- function(pars, prior_par, par_index, y_1, y_2, t, id) {

    # Order: IS, NREM, REM
    init_logit = c( 1, exp(pars[par_index$pi_logit][1]), 
                    exp(pars[par_index$pi_logit][2]), exp(pars[par_index$pi_logit][3]))

    init = init_logit / sum(init_logit)

    resp_fnc = matrix(c(1, exp(pars[par_index$misclass][1]), exp(pars[par_index$misclass][2]), 0,
                        exp(pars[par_index$misclass][3]), 1, exp(pars[par_index$misclass][4]), 0,
                        exp(pars[par_index$misclass][5]), exp(pars[par_index$misclass][6]), 1, 0,
                        0, 0, 0, 1),
                        ncol=4, byrow=TRUE)

    resp_fnc = resp_fnc / rowSums(resp_fnc)
    
    lambda_mat = matrix(c(pars[par_index$l_delta], pars[par_index$l_theta],
                          pars[par_index$l_alpha], pars[par_index$l_beta]),
                        nrow = 4)
    colnames(lambda_mat) = c("delta", "theta", "alpha", "beta")
    rownames(lambda_mat) = c("IS", "NREM", "REM", "LIMBO")

    lambda_mat = exp(lambda_mat)

    beta <- pars[par_index$beta]
  
    log_total_val = foreach(i=unique(id), .combine='+', .export = "Q", 
                            .packages = c("gtools", "expm")) %dopar% {
        
        f_i = val = 1
        y_1_i = y_1[id == i]
        y_2_i = y_2[id == i, ,drop = F] # The four proportions of waves

        t_i = t[id == i]
        disc_t_i = t_i

        d_1 = ddirichlet(x = y_2_i[1,], alpha = lambda_mat[1,])
        d_2 = ddirichlet(x = y_2_i[1,], alpha = lambda_mat[2,])
        d_3 = ddirichlet(x = y_2_i[1,], alpha = lambda_mat[3,])
        d_4 = ddirichlet(x = y_2_i[1,], alpha = lambda_mat[4,])
        
        if(y_1_i[1] <= 3) { # observed
          f_i = init %*% diag(c(d_1,d_2,d_3,0) * resp_fnc[, y_1_i[1]])
        } else { 
          f_i = init %*% diag(c(d_1,d_2,d_3,d_4))
        }

        log_norm = 0
        
        for(k in 2:length(t_i)) {
          
            P = expm((t_i[k] - t_i[k-1]) * Q(disc_t_i[k-1], beta))

            d_1 = ddirichlet(x = y_2_i[k,], alpha = lambda_mat[1,])
            d_2 = ddirichlet(x = y_2_i[k,], alpha = lambda_mat[2,])
            d_3 = ddirichlet(x = y_2_i[k,], alpha = lambda_mat[3,])
            d_4 = ddirichlet(x = y_2_i[k,], alpha = lambda_mat[4,])

            if(y_1_i[k] <= 3) { # observed
              D_i = diag(c(d_1,d_2,d_3,0) * resp_fnc[, y_1_i[k]])
            } else { # unknown 
              D_i = diag(c(d_1,d_2,d_3,d_4))
            }

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
mcmc_routine = function( y_1, y_2, t, id, init_par, prior_par, par_index,
                         steps, burnin, n_cores, ind, trialNum){

  cl <- makeCluster(n_cores, outfile="")
  registerDoParallel(cl)

  pars = init_par
  n = length(y_1)
  n_par = length(pars)
  chain = matrix( 0, steps, n_par)

  group = list(c(par_index$beta, par_index$misclass, par_index$pi_logit),
               c(par_index$l_delta, par_index$l_theta, par_index$l_alpha, 
                 par_index$l_beta))
  n_group = length(group)

  pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(group[[j]]))
  # load(paste0('Model_out/pcov_5_', trialNum - 1, '.rda'))
  pscale = rep( .0001, n_group)
  # load(paste0('Model_out/pscale_5_', trialNum - 1, '.rda'))
  accept = rep( 0, n_group)

  # Evaluate the log_post of the initial pars
  log_post_prev = fn_log_post_continuous( pars, prior_par, par_index, y_1, y_2, t, id)

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
                                        par_index, y_1, y_2, 
                                        t, id)
      print("Likelihood Evaluation:")
      print(log_post)

      # Only propose valid parameters during the burnin period
      if(ttt < burnin){
        while(!is.finite(log_post)){
          print('bad proposal')
          proposal = pars
          proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],
                                     sigma=pcov[[j]]*pscale[j])
          log_post = fn_log_post_continuous(proposal, prior_par, 
                                            par_index, y_1, y_2, 
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

  save(pcov, file = paste0("Model_out/pcov_", ind, "_", trialNum, ".rda"))
  save(pscale, file = paste0("Model_out/pscale_", ind, "_", trialNum, ".rda"))

  return(list( chain=chain[burnin:steps,], accept=accept/(steps-burnin),
               pscale=pscale))
}
# -----------------------------------------------------------------------------
