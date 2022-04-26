library(tidyverse)

fft_fnc <- function(test) {
  
  # Hann weighting window
  K = nrow(test)
  k = 1:K
  w_k = 0.5 - 0.5 * cos(2 * pi * k / (K-1))
  test$eeg = test$eeg * w_k
  # test$emg = test$emg * w_k
  
  # Formatting
  omega = seq(0, 5000, length.out =nrow(test))
  
  omega_color = rep(5, length(omega))
  
  omega_color[which(omega >= 0.8   & omega <= 4)] = 1
  omega_color[which(omega >= 4.2   & omega <= 8)] = 2
  omega_color[which(omega >= 8.2   & omega <= 13)] = 3
  omega_color[which(omega >= 13.2  & omega <= 20)] = 4
  
  fft_new = fft(test$eeg)
  # fft_new = fft(test$emg)
  
  fft_delta = fft_new[omega_color == 1]
  fft_theta = fft_new[omega_color == 2]
  fft_alpha = fft_new[omega_color == 3]
  fft_beta  = fft_new[omega_color == 4]
  
  s_power_delta = sum(abs(fft_delta)^2) 
  s_power_theta = sum(abs(fft_theta)^2) 
  s_power_alpha = sum(abs(fft_alpha)^2) 
  s_power_beta  = sum(abs(fft_beta)^2)  
  
  # return(c(s_power_delta, s_power_theta, s_power_alpha, s_power_beta))
  
  s_power_delta = sum(abs(fft_delta)^2) / (2 * pi * length(fft_delta))
  s_power_theta = sum(abs(fft_theta)^2) / (2 * pi * length(fft_theta))
  s_power_alpha = sum(abs(fft_alpha)^2) / (2 * pi * length(fft_alpha))
  s_power_beta  = sum(abs(fft_beta)^2)  / (2 * pi * length(fft_beta))
  
  s_power_total = s_power_delta + s_power_theta + s_power_alpha + s_power_beta
  
  s_pow_delta_norm = s_power_delta / s_power_total
  s_pow_theta_norm = s_power_theta / s_power_total
  s_pow_alpha_norm = s_power_alpha / s_power_total
  s_pow_beta_norm = s_power_beta / s_power_total
  
  return(c(s_pow_delta_norm, s_pow_theta_norm, s_pow_alpha_norm, s_pow_beta_norm))
}

fft_fnc_30 <- function(test) {
  
  # # test is 30 seconds of data
  # pow_30 = matrix(nrow = 6, ncol = 4)
  
  # test is 15 seconds of data
  pow_30 = matrix(nrow = 3, ncol = 4)
  
  colnames(pow_30) = c("Delta", "Theta", "Alpha", "Beta")
  
  # seq_30 = seq(0, 150000, 25000)
  seq_30 = seq(0, 75000, 25000)
  
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

mice_format = data.frame("t1"    = NULL, 
                         "t2"    = NULL,
                         "state" = NULL,
                         "delta" = NULL,
                         "theta" = NULL,
                         "alpha" = NULL,
                         "beta"  = NULL,
                         "ptnum" = NULL)

state_names <- c("limbo", "is", "nrem", "rem", "")

pat_num = 1

for (ll in 1:length(Sys.glob("Data_format/Big_data_format/*"))) {
  
  Dir <- Sys.glob("Data_format/Big_data_format/*")[ll]
  print(ll)
  print(Dir)
  
  if(ll < 8) {
    mice_data = read.csv(Dir) 
  } else {
    load(Dir)
    mice_data = mice_data_temp
  }
  
  if(length(unique(mice_data$state)) == 1) {
    if (unique(mice_data$state) == "") {
      print("No informative state")
      next
    }
  }
  
  mice_data$eeg = mice_data$eeg - mean(mice_data$eeg)
  mice_data$emg = mice_data$emg - mean(mice_data$emg)
  
  temp_df = data.frame("t1" = NA, "t2" = NA, "state" = NA, "delta" = NA,
                       "theta" = NA, "alpha" = NA, "beta" = NA)
  
  # How do these proportions change over time (every 5 seconds)
  index_seq = seq(0, nrow(mice_data), 5000)
  
  # How do these proportions change over time (every 15 seconds)
  # index_seq = seq(0, nrow(mice_data), 75000)
  
  # # How do these proportions change over time (every 30 seconds)
  # index_seq = seq(0, nrow(mice_data), 150000)
  
  for(i in 1:(length(index_seq) - 1)) {
    
    s_ind = index_seq[i] + 1
    e_ind = index_seq[i+1]
    test = mice_data[s_ind:e_ind, ] 
    
    wave_prop = fft_fnc(test)
    # wave_prop = fft_fnc_30(test)
    
    t1 = s_ind
    t2 = e_ind
    
    s = ""
    
    for(jj in 1:length(unique(test$state))) {
      temp_s = unique(test$state)[jj]
      state_num = which(state_names == temp_s)
      
      if(state_num == 5) {state_num = 99}
      s = paste0(s, state_num)
    }
    
    temp_df[i,] = c(as.numeric(t1), as.numeric(t2), 
                    as.numeric(s), as.numeric(wave_prop))
  }
  
  temp_df$ptnum = rep(pat_num, nrow(temp_df))
  temp_df$state[which(temp_df$state > 4)] = 99
  pat_num = pat_num + 1
  
  # Center and scale the time coefficient
  temp_df$t1 = temp_df$t1 - mean(temp_df$t1)
  temp_df$t1 = temp_df$t1 / sd(temp_df$t1)
  
  mice_format = rbind(mice_format, temp_df)
  
  print(head(temp_df, 3))
  print(table(temp_df$state))
}

rownames(mice_format) = NULL

save(mice_format, file = "Data_format/mice_format_total_new_5sec.rda")


# New histograms
load('Data_format/mice_format_total_new_15sec.rda')
# Looking into the 5 second distribution
temp = mice_format[mice_format$delta < 0.2, ]
t_full = as.data.frame(table(mice_format$ptnum))
t_sub = as.data.frame(table(temp$ptnum))
t_sub$Freq / t_full$Freq
which((t_sub$Freq / t_full$Freq) < 0.9) #c(1,3,4,5,6)

mice_format = mice_format[!(mice_format$ptnum %in% c(4,23)), ]
# mice_format = mice_format[mice_format$ptnum == 1, ]
# 
pdf('Plots/Supplement/mice_format_total_new_15_filt_3.pdf')
par(mfrow = c(2,2))
hist(mice_format$delta[mice_format$state == 2],
     main = "Delta, state 2")
hist(mice_format$delta[mice_format$state == 3],
     main = "Delta, state 3")
hist(mice_format$delta[mice_format$state == 4],
     main = "Delta, state 4")
hist(mice_format$delta[mice_format$state == 99],
     main = "Delta, state 99")
hist(mice_format$theta[mice_format$state == 2],
     main = "Theta, state 2")
hist(mice_format$theta[mice_format$state == 3],
     main = "Theta, state 3")
hist(mice_format$theta[mice_format$state == 4],
     main = "Theta, state 4")
hist(mice_format$theta[mice_format$state == 99],
     main = "Theta, state 99")
hist(mice_format$alpha[mice_format$state == 2],
     main = "alpha, state 2")
hist(mice_format$alpha[mice_format$state == 3],
     main = "alpha, state 3")
hist(mice_format$alpha[mice_format$state == 4],
     main = "alpha, state 4")
hist(mice_format$alpha[mice_format$state == 99],
     main = "alpha, state 99")
hist(mice_format$beta[mice_format$state == 2],
     main = "beta, state 2")
hist(mice_format$beta[mice_format$state == 3],
     main = "beta, state 3")
hist(mice_format$beta[mice_format$state == 4],
     main = "beta, state 4")
hist(mice_format$beta[mice_format$state == 99],
     main = "beta, state 99")

dev.off()

for (ll in 1:length(Sys.glob("Data_format/Exported_ephys_and_state/*.csv"))) {
  
  Dir <- Sys.glob("Data_format/Exported_ephys_and_state/*")[ll]
  mice_data = read.csv(Dir) 
  print(paste0("Mean: ", round(mean(mice_data$eeg), 8), " SD: ", sd(mice_data$eeg)))
  print("Center")
  mice_data$eeg = mice_data$eeg - mean(mice_data$eeg)
  print(paste0("Mean: ", round(mean(mice_data$eeg), 8), " SD: ", sd(mice_data$eeg)))
  print("Scaling")
  mice_data$eeg = mice_data$eeg / sd(mice_data$eeg)
  print(paste0("Mean: ", round(mean(mice_data$eeg), 8), " SD: ", sd(mice_data$eeg)))
  
  
}
