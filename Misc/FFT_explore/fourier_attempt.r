library(signal)
library(eegkit)

# Exploring the sparcity of the state labels
# for(j in 1:nrow(mice_data)) {
#   if(mice_data$state[j] != "<undefined>") {
#     print(paste0("Time: ", mice_data$t[j], ", State: ", mice_data$state[j]))
#   }
# }

Dir <- "~/Dropbox/Shared_HMM_ICU/mouse_data/"

file_names = c("Mice_Data/WT 08 20210309 03 Penetrating Arteriole 064 ECoG, EMG and sleep.csv",
               "Mice_Data/WT 08 20210311 06 Penetrating Arteriole 065 ECoG, EMG and sleep.csv",
               "Mice_Data/WT 08 20210311 10 Penetrating Arteriole 066 ECoG, EMG and sleep.csv",
               "Mice_Data/WT 09 20210308 03 Penetrating Arteriole 073 ECoG, EMG and sleep.csv")

mice_data = read.csv(file_names[2])

# We will get these estimates every 5 seconds
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

# How do these proportions change over time (every 5 seconds)
index_seq = seq(1, nrow(mice_data), 2560)
prop_change = data.frame("t1" = rep(NA, length(index_seq) - 1), 
                         "t2" = rep(NA, length(index_seq) - 1),
                         "state" = rep(NA, length(index_seq) - 1),
                         "delta" = rep(NA, length(index_seq) - 1),
                         "theta" = rep(NA, length(index_seq) - 1),
                         "alpha" = rep(NA, length(index_seq) - 1),
                         "beta" = rep(NA, length(index_seq) - 1))

for(i in 1:(length(index_seq) - 1)) {
  
  test = mice_data[index_seq[i]:(index_seq[i+1] - 1), ] # grabs data every 5 seconds
  
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
  
  prop_change[i,] = c(t1, t2, s, wave_prop)
}

plot(prop_change$t1, prop_change$delta, type = "l", lty = 1, lwd = 2, xlim = c(0,300))
lines(prop_change$t1, prop_change$theta, col = "blue", lty = 1, lwd = 2)
lines(prop_change$t1, prop_change$alpha, col = "red", lty = 1, lwd = 2)
lines(prop_change$t1, prop_change$beta, col = "green", lty = 1, lwd = 2)

plot(prop_change$t1, prop_change$delta, col = as.character(prop_change$state))


# NOTE on 'eegkit' package -------------------------
# * the wiki for "passband" has a good visual on what the filtering is doing
#
#
#
#
# NOTE on missing labels ---------------------------
# DF1: t = 287 sec is the last labeled
# DF2: t = 278 sec
# DF3: t = 289 sec
# DF4: t = 290 sec
#
