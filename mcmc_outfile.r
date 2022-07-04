requiredPackages = c('tidyverse','gridExtra')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

dir = 'Model_out/' # Change this everytime!!!! ****************

# Size of posterior sample from mcmc chains
n_post = 5000
# Step number at 3ich the adaptive tuning scheme was frozen
burnin = 5000
# Total number of steps the mcmc algorithm is computed for
steps = 30000
# Matrix row indices for the posterior sample to use for GFF computation
index_post = (steps - burnin - n_post + 1):(steps - burnin)

# par_index = list( beta=1:22, misclass=23:34, pi_logit=35:36, 
#                   l_delta = 37:40, l_theta=41:44, l_alpha=45:48, l_beta=49:52)
par_index = list( beta=1:24, pi_logit=25:27, l_delta = 28:31, l_theta=32:35, 
                  l_alpha=36:39, l_beta=40:43)

index_seeds = c(1,3:5,8)
trialNum = 44 # Change this everytime!!!! ****************

# Initial parameters for the 30-s epochs
true_par = c(c(matrix(c(-6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0), ncol=2, byrow = T)),
            c(0, 0, 0),
            log(c(0.42844057, 0.52802631, 0.46033976, 0.52947184,
              0.34666468, 0.30778351, 0.36753584, 0.30819960,
              0.17758807, 0.12382649, 0.13143568, 0.11718367,
              0.04730668, 0.04036369, 0.04068871, 0.04514489)))

labels <- c("Baseline: IS --> NREM", "Baseline: IS --> REM", "Baseline: IS --> LIMBO",   
            "Baseline: NREM --> IS", "Baseline: NREM --> REM", "Baseline: NREM --> LIMBO",
            "Baseline: REM --> IS", "Baseline: REM --> NREM", "Baseline: REM --> LIMBO",
            "Baseline: LIMBO --> IS", "Baseline: LIMBO --> NREM", "Baseline: LIMBO --> REM",
            "Time: IS --> NREM", "Time: IS --> REM", "Time: IS --> LIMBO",      
            "Time: NREM --> IS",   "Time: NREM --> REM", "Time: NREM --> LIMBO",
            "Time: REM --> IS", "Time: REM --> NREM", "Time: REM --> LIMBO",
            "Time: LIMBO --> IS", "Time: LIMBO --> NREM", "Time: LIMBO --> REM",
            # "P( obs. IS | true NREM )", "P( obs. IS | true REM )",
            # "P( obs. NREM | true IS )", "P( obs. NREM | true REM )", 
            # "P( obs. REM | true IS )", "P( obs. REM | true NREM )", 
            "P( init NREM )", "P( init REM )", "P( init LIMBO )",
            "Delta (IS)", "Delta (NREM)", "Delta (REM)", "Delta (LIMBO)",
            "Theta (IS)", "Theta (NREM)", "Theta (REM)", "Theta (LIMBO)",
            "Alpha (IS)", "Alpha (NREM)", "Alpha (REM)", "Alpha (LIMBO)",
            "Beta (IS)", "Beta (NREM)", "Beta (REM)", "Beta (LIMBO)")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Calculating Credible Sets ---------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
cred_set = vector(mode = 'list', length = length(true_par))
for(i in 1:length(cred_set)) { cred_set[[i]] = data.frame('lower' = c(-1), 'upper' = c(-1)) }

ind = 0

for (i in index_seeds) {
    file_name = paste0(dir,'mcmc_out_',toString(i),'_', trialNum,'.rda')
    if(file.exists(file_name)) {
        load(file_name)
        ind = ind + 1
        # Inverse logit to convert back to probabilities
        
        for(j in 1:length(true_par)) {
    
            cred_set[[j]][ind,1] =  round(quantile( mcmc_out$chain[index_post,j],
                                        prob=.025), 4)
            cred_set[[j]][ind,2] =  round(quantile( mcmc_out$chain[index_post,j],
                                        prob=.975), 4)
        }
    }
}

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
    # print(paste0("Coverage for parameter ", val, " is: ", covrg))
}


# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
chain_list = vector(mode = "list", length = ind)
post_means = matrix(nrow = ind, ncol = length(labels))

ind = 0

for(seed in index_seeds){
    file_name = paste0(dir,'mcmc_out_',toString(seed),'_', trialNum, '.rda')
    if (file.exists(file_name)) {
        load(file_name)
        ind = ind + 1
        print(paste0(ind, ": ", file_name))
        print(mcmc_out$accept)
        in_post_temp = tail(index_post, 300)
        print(colMeans(mcmc_out$chain[in_post_temp,]))
        
        chain_list[[ind]] = mcmc_out$chain[index_post,]
        # post_means[ind,] <- colMeans(mcmc_out$chain[index_post,])
    }
}

# print(post_means)

# Plot and save the mcmc trace plots and histograms.
pdf(paste0('Plots/Seed_5_Progression/mcmc_total', '_', trialNum, '.pdf')) # *****
par(mfrow=c(4, 2))

stacked_chains = do.call( rbind, chain_list)
par_mean = par_median = upper = lower = rep( NA, length(labels))
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
    # abline( v=true_par[r], col='green', lwd=2, lty=2)
    abline( v=lower[r], col='purple', lwd=2, lty=2)

}


# for(r in 1:length(labels)) {
#     # Adding the boxplots
#     yVar = disc_type = x_label = NULL
    
#     yVar = post_means[[1]][,r]
#     disc_type = rep("Continuous", nrow(post_means[[1]]))
#     x_label = paste0("Coverage is: ", round(cov_df[r], digits=3))

#     plot_df = data.frame(yVar = yVar, disc_type = disc_type)
#     VP[[r]] = ggplot(plot_df, aes(x=disc_type, y = yVar)) +
#       geom_violin(trim=FALSE) +
#       geom_boxplot(width=0.1) +
#       ggtitle(labels[r]) +
#       ylab(paste0("Parameter Value: ", round(true_par[r], 3))) +
#       xlab(x_label) +
#       geom_hline(yintercept=true_par[r], linetype="dashed", color = "red") +
#       theme(text = element_text(size = 7))
# }

# grid.arrange(VP[[1]], VP[[2]], VP[[3]], VP[[4]], VP[[5]],
#              VP[[6]], VP[[7]], VP[[8]], VP[[9]], ncol=3, nrow =3)
# grid.arrange(VP[[10]], VP[[11]], VP[[12]], VP[[13]], VP[[14]],
#              VP[[15]], VP[[16]], VP[[17]], VP[[18]], ncol=3, nrow =3)
# grid.arrange(VP[[19]], VP[[20]], VP[[21]], ncol=3, nrow =3)

dev.off()

# save(post_means, file = paste0("Plots/post_means_", ind_num, ".rda"))
# save(cov_df, file = paste0("Plots/cov_df_", ind_num, ".rda"))
