# args = commandArgs(TRUE)

fft_fnc <- function(test) {
  
  # Hann weighting window
  K = nrow(test)
  k = 1:K
  w_k = 0.5 - 0.5 * cos(2 * pi * k / (K-1))
  test$ecog = test$ecog * w_k
  # test$emg = test$emg * w_k

  # Formatting
  omega = seq(0, 1/(test$t[2] - test$t[1]), length.out =nrow(test))
  
  omega_color = rep(5, length(omega))
  
  omega_color[which(omega >= 0.8 & omega <= 4 )] = 1
  omega_color[which(omega >= 4.2   & omega <= 8)] = 2
  omega_color[which(omega >= 8.2   & omega <= 13)] = 3
  omega_color[which(omega >= 13.2  & omega <= 20)] = 4
  
  fft_new = fft(test$ecog)
  # fft_new = fft(test$emg)
  
  fft_delta = fft_new[omega_color == 1]
  fft_theta = fft_new[omega_color == 2]
  fft_alpha = fft_new[omega_color == 3]
  fft_beta  = fft_new[omega_color == 4]
  
  s_power_delta = sum(abs(fft_delta)^2) 
  s_power_theta = sum(abs(fft_theta)^2) 
  s_power_alpha = sum(abs(fft_alpha)^2) 
  s_power_beta  = sum(abs(fft_beta)^2)  
  
  return(c(s_power_delta, s_power_theta, s_power_alpha, s_power_beta))

  # s_power_delta = sum(abs(fft_delta)^2) / (2 * pi * length(fft_delta))
  # s_power_theta = sum(abs(fft_theta)^2) / (2 * pi * length(fft_theta))
  # s_power_alpha = sum(abs(fft_alpha)^2) / (2 * pi * length(fft_alpha))
  # s_power_beta  = sum(abs(fft_beta)^2)  / (2 * pi * length(fft_beta))
  
  # s_power_total = s_power_delta + s_power_theta + s_power_alpha + s_power_beta
  
  # s_pow_delta_norm = s_power_delta / s_power_total; s_pow_delta_norm
  # s_pow_theta_norm = s_power_theta / s_power_total; s_pow_theta_norm
  # s_pow_alpha_norm = s_power_alpha / s_power_total; s_pow_alpha_norm
  # s_pow_beta_norm = s_power_beta / s_power_total; s_pow_beta_norm
  
  # return(c(s_pow_delta_norm, s_pow_theta_norm, s_pow_alpha_norm, s_pow_beta_norm))
}

fft_fnc_30 <- function(test) {
  
  # test is 30 seconds of data
  pow_30 = matrix(nrow = 6, ncol = 4)
  colnames(pow_30) = c("Delta", "Theta", "Alpha", "Beta")

  seq_30 = seq(0, 15360, 2560)

  for(i in 2:length(seq_30)) {
    start_row = seq_30[i-1] + 1
    end_row = seq_30[i]
    temp_5 = test[start_row:end_row, ]

    pow_30[i-1, ] = fft_fnc(temp_5)
  }
  
  power_band = colMeans(pow_30)

  s_power_total = sum(power_band)
  
  s_pow_delta_norm = power_band[1] / s_power_total
  s_pow_theta_norm = power_band[2] / s_power_total
  s_pow_alpha_norm = power_band[3] / s_power_total
  s_pow_beta_norm =  power_band[4] / s_power_total
  
  return(c(s_pow_delta_norm, s_pow_theta_norm, s_pow_alpha_norm, s_pow_beta_norm))
}

Dir <- "~/Dropbox/Shared_HMM_ICU/mouse_data/WT 08 20210309 03 Penetrating Arteriole 064 ECoG, EMG and sleep.csv"

mice_data = read.csv(Dir)

state_names <- c("Limbo", "Clean IS", "Clean NREM", "Clean REM", "<undefined>")

# # How do these proportions change over time (every 5 seconds)
# index_seq = seq(0, nrow(mice_data), 2560)

# How do these proportions change over time (every 30 seconds)
index_seq = seq(0, nrow(mice_data), 15360)

mice_format = data.frame("t1" = rep(NA, length(index_seq) - 1), 
                         "t2" = rep(NA, length(index_seq) - 1),
                         "state" = rep(NA, length(index_seq) - 1),
                         "delta" = rep(NA, length(index_seq) - 1),
                         "theta" = rep(NA, length(index_seq) - 1),
                         "alpha" = rep(NA, length(index_seq) - 1),
                         "beta" = rep(NA, length(index_seq) - 1))

for(i in 1:(length(index_seq) - 1)) {
  
  s_ind = index_seq[i] + 1
  e_ind = index_seq[i+1]
  test = mice_data[s_ind:e_ind, ] # grabs data every 30 seconds
  
  wave_prop = fft_fnc_30(test)
  # wave_prop = fft_fnc(test)
  
  t1 = head(test$t,1)
  t2 = tail(test$t,1)
  
  s = ""
  
  for(jj in 1:length(unique(test$state))) {
    temp_s = unique(test$state)[jj]
    state_num = which(state_names == temp_s)
    
    if(state_num == 5) {state_num = 99}
    s = paste0(s, state_num)
  }
  
  print(paste0(i, "  ", s))
  mice_format[i,] = c(as.numeric(t1), as.numeric(t2), 
                      as.numeric(s), as.numeric(wave_prop))
}

mice_format$ptnum = rep(1, nrow(mice_format))
mice_format$state[which(mice_format$state > 4)] = 99

# # Center and scale the time coefficient
mice_format$t1 = mice_format$t1 - mean(mice_format$t1)
mice_format$t1 = mice_format$t1 / sd(mice_format$t1)

print(unique(mice_format$state))
save(mice_format, file = "Data_format/mice_format.rda")

