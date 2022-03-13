requiredPackages = c('tidyverse','gridExtra')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

nFrames = 1

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
            "P( init IS )", "P( init NREM )")

mice_data = matrix(data=-1, nrow = nFrames, ncol = length(labels))

for (i in 1:nFrames) {

  load(paste0("Model_out/Output_msm", i, ".rda"))
  mice_data[i,] = Output_msm$opt$par

}

meanValues <- colMeans(mice_data)


pdf("Plots/msm.pdf", onefile = T)
VP <- vector(mode="list", length = length(labels))
for(i in 1:length(labels)) {
    plot_df = data.frame(yVar = mice_data[,i])
    VP[[i]] = ggplot(plot_df, aes(x="", y = yVar)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.1) +
      ggtitle(labels[i]) +
      ylab(paste0("Parameter Mean: ", meanValues[i])) +
      xlab('') +
      geom_hline(yintercept=meanValues[i], linetype="dashed", color = "red") +
      theme(text = element_text(size = 7))

}
grid.arrange(VP[[1]], VP[[2]], VP[[3]], VP[[4]], VP[[5]],
             VP[[6]], VP[[7]], VP[[8]], VP[[9]], ncol=3, nrow =3)
grid.arrange(VP[[10]], VP[[11]], VP[[12]], VP[[13]],
             VP[[14]], VP[[15]], VP[[16]], ncol=3, nrow =3)
dev.off()
