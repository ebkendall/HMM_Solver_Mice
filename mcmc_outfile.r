requiredPackages = c('tidyverse','gridExtra')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

# args = commandArgs(TRUE)
# ind_num = as.numeric(args[1])
# ind_t = as.numeric(args[2])
ind_num = 1
ind_t = 30 # or 5

dir = 'Model_out/'

# Size of posterior sample from mcmc chains
n_post = 5000
# Step number at 3ich the adaptive tuning scheme was frozen
burnin = 5000
# Total number of steps the mcmc algorithm is computed for
steps = 10000
# Matrix row indices for the posterior sample to use for GFF computation
index_post = (steps - burnin - n_post + 1):(steps - burnin)

par_index = list( beta=1:22, misclass=23:34, pi_logit=35:36, 
                  l_delta = 37:40, l_theta=41:44, l_alpha=45:48, l_beta=49:52)

index_seeds = 1:ind_num


# Initial parameters for the 5-s epochs
# true_par = c(c(matrix(c(2.071653   , 2.85576,
#                         1.212561   , -3.418799,
#                         -4.198075  , 0.5761976,
#                         -1.882976  , 4.516193,
#                         -0.01241128, -2.303688,
#                         -0.6022949 , 3.432778,
#                         -0.9710429 , -2.282282,
#                         1.126578   , 4.04868,
#                         -0.6175124 , 0.6705703,
#                         3.312069   , 0.6658912,
#                         -2.792293  , 0.1919586), ncol=2, byrow=T)),
#                         c(-9.266642, 0.9177677, -4.659014,
#                           -4.381079, 2.809255, -2.815062,
#                           -1.498909, -2.702002, -2.519476,
#                           -5.202432, -2.634818, -1.023371),
#                         c(-2.573828,  1.353468), 
#                         c(2.406984, 7.39707, 0.006235556, 0.9011201),
#                         c(2.526422, 2.422841, 1.793587, 0.3687743),
#                         c(-1.534735, 0.01483018, 1.649471, 5.28734),
#                         c(3.082263, 2.136501, 1.184794, -0.081426944))

# Initial parameters for the 30-s epochs
true_par = c(c(matrix(c(4.101165   , 2.936901,
                        2.592022   , 0.4181154,
                        -3.341068  , -1.838769,
                        -2.744266  , 6.732399,
                        1.360784   , -2.182402,
                        -0.4102578 , 1.793358,
                        -2.385557  , 1.524106,
                        1.001188   , 4.97785,
                        -4.059381  , -1.370959,
                        -0.4864039 , 3.183361,
                        -2.636722  , 1.41935), ncol=2, byrow=T)),
                        c(-7.736442, 4.381084, -2.01876,
                          -1.401513, 2.879382, -4.285909,
                          -1.589179, -0.1894496, -4.367929,
                          -7.096464, -4.536896, -2.845439),
                        c(-1.049866, 1.229301), 
                        c(3.737221, 5.801327, -0.6352715, -2.009164),
                        c(3.5487, 3.786715, 2.98324, 2.350825),
                        c(-3.159109, 2.168725, 4.431401, 4.442483),
                        c(3.644593, 3.197575, 2.481969, 1.792291))

labels <- c("Baseline: LIMBO --> IS",  "Baseline: LIMBO --> NREM",
            "Baseline: IS --> LIMBO",  "Baseline: IS --> NREM",
            "Baseline: IS --> REM",    "Baseline: NREM --> LIMBO",
            "Baseline: NREM --> IS",   "Baseline: NREM --> REM",
            "Baseline: REM --> LIMBO", "Baseline: REM --> IS",
            "Baseline: REM --> NREM", 
            "Time: LIMBO --> IS",  "Time: LIMBO --> NREM",
            "Time: IS --> LIMBO",  "Time: IS --> NREM",
            "Time: IS --> REM",    "Time: NREM --> LIMBO",
            "Time: NREM --> IS",   "Time: NREM --> REM",
            "Time: REM --> LIMBO", "Time: REM --> IS",
            "Time: REM --> NREM",
            "P( obs. LIMBO | true IS )", "P( obs. LIMBO | true NREM )",
            "P( obs. LIMBO | true REM )", "P( obs. IS | true LIMBO )", 
            "P( obs. IS | true NREM )", "P( obs. IS | true REM )",
            "P( obs. NREM | true LIMBO )", "P( obs. NREM | true IS )",
            "P( obs. NREM | true REM )", "P( obs. REM | true LIMBO )",
            "P( obs. REM | true IS )", "P( obs. REM | true NREM )",
            "P( init IS )", "P( init NREM )", 
            "Delta (LIMBO)", "Theta (LIMBO)", "Alpha (LIMBO)", "Beta (LIMBO)",
            "Delta (IS)", "Theta (IS)", "Alpha (IS)", "Beta (IS)",
            "Delta (NREM)", "Theta (NREM)", "Alpha (NREM)", "Beta (NREM)",
            "Delta (REM)", "Theta (REM)", "Alpha (REM)", "Beta (REM)")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Calculating Credible Sets ---------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
cred_set = vector(mode = 'list', length = length(true_par))
for(i in 1:length(cred_set)) { cred_set[[i]] = data.frame('lower' = c(-1), 'upper' = c(-1)) }

ind = 0

for (i in index_seeds) {
    file_name = paste0(dir,'mcmc_out_',toString(i),'_10_', ind_t,'.rda')

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
    print(paste0("Coverage for parameter ", val, " is: ", covrg))
}


# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------
chain_list = vector(mode = "list", length = ind)
post_means = matrix(nrow = ind, ncol = length(labels))

ind = 0

for(seed in index_seeds){
    file_name = paste0(dir,'mcmc_out_',toString(seed),'_10_', ind_t, '.rda')
    if (file.exists(file_name)) {
        load(file_name)
        ind = ind + 1
        print(mcmc_out$accept)

        chain_list[[ind]] = mcmc_out$chain[index_post,]
        post_means[ind,] <- colMeans(mcmc_out$chain[index_post,])
    }
}

print(post_means)

# Plot and save the mcmc trace plots and histograms.
pdf(paste0('Plots/mcmc_', ind_num, '_10_', ind_t, '.pdf'))
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
    abline( v=true_par[r], col='green', lwd=2, lty=2)
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
