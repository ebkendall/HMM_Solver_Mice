requiredPackages = c('tidyverse','gridExtra')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

# args = commandArgs(TRUE)
# ind_num = as.numeric(args[1])
# ind_t = as.numeric(args[2])
ind_num = 50

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


# Initial parameters for the 30-s epochs
true_par = c(c(matrix(c(-2.26568339,  0.08766060,
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
    file_name = paste0(dir,'mcmc_out_',toString(i),'_', steps/1000, '.rda')

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
    file_name = paste0(dir,'mcmc_out_',toString(seed),'_', steps/1000, '.rda')
    if (file.exists(file_name)) {
        load(file_name)
        ind = ind + 1
        print(mcmc_out$accept)

        chain_list[[ind]] = mcmc_out$chain[index_post,]
        post_means[ind,] <- colMeans(mcmc_out$chain[index_post,])
    }
}

# print(post_means)

# Plot and save the mcmc trace plots and histograms.
pdf(paste0('Plots/mcmc_sim', '_', steps/1000, '.pdf'))
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
