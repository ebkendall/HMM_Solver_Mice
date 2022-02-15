# This is where we will load the time series and clean the data such that it 
# is formatted in the way that we want. In other words, for every 5 seconds,
# we will have a time measure, and then four proportions for the different times
# as well as the different states.


# Steps
# 1) load the existing data
# 2) split the data into the respective time intervals
# 3) run the eeg splitting technique such that we have 4 proportions
# 4) the dataframe should look like
#         colnames(df) = c("t", "delta_prop", "theta_prop",
#                          "alpha_prop", "beta_prop", "state",
#                          ... other covariates)

# Questions / Concerns
# 1) How to deal with the fact that separate data frames are involving separate
#   states (i.e. one has sleep properties and the other has awake properties)
# 2) How to deal with time periods with overlapping states

# Good example of SLEEP: 
#    WT 08 20210309 03 Penetrating Arteriole 064 ECoG, EMG and sleep.csv
# Good example of AWAKE: 
#    WT 11 20210705 05 Penetrating Arteriole 114 ECoG, EMG and sleep.csv

library(signal)
library(eegkit)

fft_fnc <- function(test) {
  
  f_delta = eegfilter(test$ecog, Fs = length(which(mice_data$t < 1)), lower = 1, upper = 4, method = "butter")
  f_theta = eegfilter(test$ecog, Fs = length(which(mice_data$t < 1)), lower = 4, upper = 7, method = "butter")
  f_alpha = eegfilter(test$ecog, Fs = length(which(mice_data$t < 1)), lower = 7, upper = 12, method = "butter")
  f_beta  = eegfilter(test$ecog, Fs = length(which(mice_data$t < 1)), lower = 12, upper = 30, method = "butter")
  
  # plot(test$t, test$ecog, type = "l", lty = 1, lwd = 2, ylim = c(-1, 1))
  # lines(test$t, f_delta, col = "blue", lty = 2, lwd = 2)
  # lines(test$t, f_theta, col = "red", lty = 2, lwd = 2)
  # lines(test$t, f_alpha, col = "green", lty = 2, lwd = 2)
  # lines(test$t, f_beta, col = "purple", lty = 2, lwd = 2)
  
  
  fft_delta = fft(f_delta)
  fft_theta = fft(f_theta)
  fft_alpha = fft(f_alpha)
  fft_beta  = fft(f_beta)
  
  s_power_delta = sum(abs(fft_delta)^2) / (2 * pi * length(fft_delta))
  s_power_theta = sum(abs(fft_theta)^2) / (2 * pi * length(fft_theta))
  s_power_alpha = sum(abs(fft_alpha)^2) / (2 * pi * length(fft_alpha))
  s_power_beta  = sum(abs(fft_beta)^2)  / (2 * pi * length(fft_beta))
  
  s_power_total = s_power_delta + s_power_theta + s_power_alpha + s_power_beta
  
  s_pow_delta_norm = s_power_delta / s_power_total; s_pow_delta_norm
  s_pow_theta_norm = s_power_theta / s_power_total; s_pow_theta_norm
  s_pow_alpha_norm = s_power_alpha / s_power_total; s_pow_alpha_norm
  s_pow_beta_norm = s_power_beta / s_power_total; s_pow_beta_norm
  
  return(c(s_pow_delta_norm, s_pow_theta_norm, s_pow_alpha_norm, s_pow_beta_norm))
}

Dir <- "~/Dropbox/Shared_HMM_ICU/mouse_data/"

file_names <- c("WT 08 20210309 03 Penetrating Arteriole 064 ECoG, EMG and sleep.csv",
                "WT 11 20210705 05 Penetrating Arteriole 114 ECoG, EMG and sleep.csv")

seed_num = 1 # or 2
set.seed(seed_num)

mice_data = read.csv(paste0(Dir, file_names[seed_num]))

index_seq = seq(1, nrow(mice_data), 2560)
mice_format = data.frame("t1" = rep(NA, length(index_seq) - 1), 
                         "t2" = rep(NA, length(index_seq) - 1),
                         "state" = rep(NA, length(index_seq) - 1),
                         "delta" = rep(NA, length(index_seq) - 1),
                         "theta" = rep(NA, length(index_seq) - 1),
                         "alpha" = rep(NA, length(index_seq) - 1),
                         "beta" = rep(NA, length(index_seq) - 1))

for(i in 1:(length(index_seq) - 1)) {
  
  test = mice_data[index_seq[i]:index_seq[i+1], ] # grabs data every 5 seconds
  
  wave_prop = fft_fnc(test)
  
  t1 = head(test$t,1)
  t2 = tail(test$t,1)
  
  s = NULL
  
  print(paste0(i, "  ", unique(test$state)))
  
  if(length(unique(test$state)) == 1) {
    if(unique(test$state) == "<undefined>") {
      s = -1
    } else { s = unique(test$state) }
  } else {
    s = -1 * length(unique(test$state))
  }
  
  mice_format[i,] = c(t1, t2, s, wave_prop)
}

plot(mice_format$t1, mice_format$delta, type = "l", lty = 1, lwd = 2, xlim = c(0,300))
lines(mice_format$t1, mice_format$theta, col = "blue", lty = 1, lwd = 2)
lines(mice_format$t1, mice_format$alpha, col = "red", lty = 1, lwd = 2)
lines(mice_format$t1, mice_format$beta, col = "green", lty = 1, lwd = 2)

plot(mice_format$t1, mice_format$delta, col = as.character(mice_format$state))

# Investigating how frequently the time between states changes
freq_df = data.frame("start" = c(1), "end" = c(1), "state" = c(1))
freq_df[1,] = c(mice_data$t[1], mice_data$t[1], mice_data$state[1])
row_num = 1
for (i in 2:nrow(mice_data)) {
  if(mice_data$state[i] != mice_data$state[i-1]) {
    freq_df$end[row_num] = mice_data$t[i-1]
    row_num = row_num + 1
    freq_df[row_num,] = c(mice_data$t[i], mice_data$t[i], mice_data$state[i])
  } else {
    freq_df$end[row_num] = mice_data$t[i]
  }
}

freq_df$start = as.numeric(freq_df$start)
freq_df$end = as.numeric(freq_df$end)


