requiredPackages = c('tidyverse','gridExtra')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

dir = paste0('Model_out/')

# Size of posterior sample from mcmc chains
n_post = 5000
# Step number at 3ich the adaptive tuning scheme was frozen
burnin = 5000
# Total number of steps the mcmc algorithm is computed for
steps = 10000
# Matrix row indices for the posterior sample to use for GFF computation
index_post = (steps - burnin - n_post + 1):(steps - burnin)

index_seeds = 1:100

true_par <- c(c(matrix(c(-2.29709805,  0.09266760, -0.56262135,
                           -1.17308794, -5.10636947, -0.96162312,
                           -1.71474254, -0.04338819,  0.83882558,
                           -2.08300714,  0.03824367, -2.75345311,
                           -2.42208380,  0.11315485,  1.76897841,
                           -1.9       ,  0.15      ,  1.1  ), ncol=3, byrow=T)),
                           c(  -5.60251814, -0.84455697, -2.56906519, -2.12629033),
                           c( -6.95125291, -7.07504453), # these may be known with certainty
                           c(10, 20, 30), #only three states
                           1) # needs to be just 1 sigma


par_index = list( beta=1:18, misclass=19:22, pi_logit=23:24, mu = 25:27, sigma = 28)
# Doing the inverse logit for true_par
# true_par[par_index$pi_logit] =
#     exp(true_par[par_index$pi_logit])/(1 + exp(true_par[par_index$pi_logit]))
# true_par[par_index$misclass[1]] =
#     exp(true_par[par_index$misclass[1]])/(1 + exp(true_par[par_index$misclass[1]]))
# true_par[par_index$misclass[2:3]] =
#     exp(true_par[par_index$misclass[2:3]])/sum(c(1, exp(true_par[par_index$misclass[2:3]])))
# true_par[par_index$misclass[4]] =
#     exp(true_par[par_index$misclass[4]])/(1 + exp(true_par[par_index$misclass[4]]))

print(true_par)

labels <- c('b.l. S1  --->   S2 ',
            'b.l. S1  --->   S3 ',
            'b.l. S2  --->   S1 ',
            'b.l. S2  --->   S3 ',
            'b.l. S3    --->   S1 ',
            'b.l. S3    --->   S2 ',
            'time S1  --->   S2 ',
            'time S1  --->   S3 ',
            'time S2  --->   S1 ',
            'time S2  --->   S3 ',
            'time S3    --->   S1 ',
            'time S3    --->   S2 ',
            'sex S1   --->   S2',
            'sex S1   --->   S3',
            'sex S2  --->   S1',
            'sex S2   --->   S3',
            'sex S3    --->   S1',
            'sex S3    --->   S2',
            'P( obs. S2 | true S1 )',
            'P( obs. S1 | true S2 )',
            'P( obs. S3 | true S2 )',
            'P( obs. S2 | true S3 )',
            'P( init S2 )','P( init S3 )',
            'Mean 1', 'Mean 2', 'Mean 3',
            'Standard Deviation')

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Calculating Credible Sets ---------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
cred_set = vector(mode = 'list', length = length(true_par))
for(i in 1:length(cred_set)) {cred_set[[i]] = data.frame('lower' = c(-1), 'upper' = c(-1))}

for (i in index_seeds) {
    file_name = paste0(dir,'mcmc_out_',toString(i),'.rda')

    load(file_name)

    # Inverse logit to convert back to probabilities
    # mcmc_out$chain[,par_index$pi_logit] =
    #     exp(mcmc_out$chain[,par_index$pi_logit])/(1 + exp(mcmc_out$chain[,par_index$pi_logit]))
    # mcmc_out$chain[,par_index$misclass[1]] =
    #     exp(mcmc_out$chain[,par_index$misclass[1]])/(1 + exp(mcmc_out$chain[,par_index$misclass[1]]))
    # mcmc_out$chain[,par_index$misclass[2:3]] = exp(mcmc_out$chain[,par_index$misclass[2:3]])/(1 +
    #                                            exp(mcmc_out$chain[,par_index$misclass[2]]) +
    #                                            exp(mcmc_out$chain[,par_index$misclass[3]]))
    # mcmc_out$chain[,par_index$misclass[4]] =
    #     exp(mcmc_out$chain[,par_index$misclass[4]])/(1 + exp(mcmc_out$chain[,par_index$misclass[4]]))

    for(j in 1:length(true_par)) {
        cred_set[[j]][i,1] =  round(quantile( mcmc_out$chain[index_post,j],
                                    prob=.025), 4)
        cred_set[[j]][i,2] =  round(quantile( mcmc_out$chain[index_post,j],
                                    prob=.975), 4)
    }
}

# save(cred_set, file = paste0('Plots/cred_set_', model_name[folder], '.rda'))

# -----------------------------------------------------------------------------
# Calculating Coverage --------------------------------------------------------
# -----------------------------------------------------------------------------
cov_df <- c()
for(i in 1:length(true_par)) {
    val = true_par[i]
    top = length(which(cred_set[[i]]$lower <= val & val <= cred_set[[i]]$upper))
    bot = nrow(cred_set[[i]])
    covrg = top/bot
    cov_df[i] = covrg
    print(paste0("Coverage for parameter ", val, " is: ", covrg))
}
# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------

post_means = matrix(nrow = length(index_seeds), ncol = length(labels))
chain_list = vector(mode = "list", length = length(index_seeds))

for(seed in index_seeds){
    file_name = paste0(dir,'mcmc_out_',toString(seed),'.rda')

    load(file_name)

    print(mcmc_out$accept)

    # Inverse logit to convert back to probabilities
    # mcmc_out$chain[,par_index$pi_logit] =
    #     exp(mcmc_out$chain[,par_index$pi_logit])/(1 + exp(mcmc_out$chain[,par_index$pi_logit]))
    # mcmc_out$chain[,par_index$misclass[1]] =
    #     exp(mcmc_out$chain[,par_index$misclass[1]])/(1 + exp(mcmc_out$chain[,par_index$misclass[1]]))
    # mcmc_out$chain[,par_index$misclass[2:3]] = exp(mcmc_out$chain[,par_index$misclass[2:3]])/(1 +
    #                                            exp(mcmc_out$chain[,par_index$misclass[2]]) +
    #                                            exp(mcmc_out$chain[,par_index$misclass[3]]))
    # mcmc_out$chain[,par_index$misclass[4]] =
    #     exp(mcmc_out$chain[,par_index$misclass[4]])/(1 + exp(mcmc_out$chain[,par_index$misclass[4]]))


    chain_list[[seed]] = mcmc_out$chain[index_post,]
    post_means[seed,] <- colMeans(mcmc_out$chain[index_post,])
}

# Plot and save the mcmc trace plots and histograms.
pdf(paste0('Plots/mcmc_cont_resp.pdf'))
par(mfrow = c(4,2))

stacked_chains = do.call( rbind, chain_list)
par_mean = par_median = upper = lower = rep( NA, ncol(stacked_chains))
VP <- vector(mode="list", length = length(labels))

for(r in 1:length(labels)){

	plot( NULL, xlab=NA, ylab=NA, main=labels[r], xlim=c(1,length(index_post)),
		    ylim=range(stacked_chains[,r]) )

	for(seed in 1:length(chain_list)) lines( chain_list[[seed]][,r], type='l', col=seed)

	par_mean[r] = round( mean(stacked_chains[,r]), 4)
	par_median[r] = round( median(stacked_chains[,r]), 4)
	upper[r] = round( quantile( stacked_chains[,r], prob=.975), 4)
	lower[r] = round( quantile( stacked_chains[,r], prob=.025), 4)

	hist( stacked_chains[,r], breaks=sqrt(nrow(stacked_chains)), ylab=NA, main=NA,
	      freq=F, xlab=paste0('Mean = ',toString(par_mean[r]),
	                           ' Median = ',toString(par_median[r])))
	abline( v=upper[r], col='red', lwd=2, lty=2)
	abline( v=true_par[r], col='green', lwd=2, lty=2)
	abline( v=lower[r], col='purple', lwd=2, lty=2)

}



for(r in 1:length(labels)) {
    # Adding the boxplots
    yVar = disc_type = x_label = NULL

    yVar = post_means[,r]
    disc_type = rep("Continuous", nrow(post_means))
    x_label = paste0("Coverage is: ", cov_df[r])

    plot_df = data.frame(yVar = yVar, disc_type = disc_type)
    VP[[r]] = ggplot(plot_df, aes(x=disc_type, y = yVar)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.1) +
      ggtitle(labels[r]) +
      ylab(paste0("Parameter Value: ", round(true_par[r], 3))) +
      xlab(x_label) +
      geom_hline(yintercept=true_par[r], linetype="dashed", color = "red") +
      theme(text = element_text(size = 7))
}

grid.arrange(VP[[1]], VP[[2]], VP[[3]], VP[[4]], VP[[5]],
             VP[[6]], VP[[7]], VP[[8]], VP[[9]], ncol=3, nrow =3)
grid.arrange(VP[[10]], VP[[11]], VP[[12]], VP[[13]], VP[[14]],
             VP[[15]], VP[[16]], VP[[17]], VP[[18]], ncol=3, nrow =3)
grid.arrange(VP[[19]], VP[[20]], VP[[21]],
             VP[[22]], VP[[23]], VP[[24]],
             VP[[25]], VP[[26]], VP[[27]], ncol=3, nrow =3)
grid.arrange(VP[[28]], ncol=3, nrow =3)

dev.off()
