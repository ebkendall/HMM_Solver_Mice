Dir <- "~/Dropbox/Shared_HMM_ICU/mouse_data/WT 08 20210309 03 Penetrating Arteriole 064 ECoG, EMG and sleep.csv"

mice_data = read.csv(Dir)
state_names <- c("Limbo", "Clean IS", "Clean NREM", "Clean REM", "<undefined>")

# Seeing if there exists a pattern for different states in the fourier series
fft_fnc <- function(test) {
  
  omega = seq(0, 1/(test$t[2] - test$t[1]), length.out =length(test$t))
  
  omega_color = omega
  
  omega_color[which(omega >= 0.8 & omega < 4 )] = 1
  omega_color[which(omega >= 4   & omega < 7 )] = 2
  omega_color[which(omega >= 7   & omega < 12)] = 3
  omega_color[which(omega >= 12  & omega < 30)] = 4
  omega_color[which(omega >= 30 | omega < 0.8)] = 5
  
  fft_new = fft(test$ecog)
  
  fft_delta = fft_new[omega_color == 1]
  fft_theta = fft_new[omega_color == 2]
  fft_alpha = fft_new[omega_color == 3]
  fft_beta  = fft_new[omega_color == 4]
  
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
  
  s = ""
  
  for(jj in 1:length(unique(test$state))) {
    temp_s = unique(test$state)[jj]
    state_num = which(state_names == temp_s)
    
    if(state_num == 5) {state_num = 99}
    s = paste0(s, state_num)
  }
  
  # if(length(unique(test$state)) == 1) {
  #   if(unique(test$state) == "<undefined>") {
  #     s = -1
  #   } else { s = unique(test$state) }
  # } else {
  #   s = -1 * length(unique(test$state))
  # }
  
  print(paste0(i, "  ", s))
  mice_format[i,] = c(as.numeric(t1), as.numeric(t2), 
                      as.numeric(s), as.numeric(wave_prop))
}

mice_format$ptnum = rep(1, nrow(mice_format))
mice_format$state[which(mice_format$state > 4)] = 7

band_type = c("Delta", "Theta", "Alpha", "Beta")

# pdf("stateChangeBand.pdf")
# par(mfrow = c(4,1))
# state_ind = 2
# for(i in sort(unique(mice_format$state))) {
#   print(i)
#   
#   color_type = mice_format$state
#   color_type[which(mice_format$state != i)] = 0
#   
#   plot(mice_format$t1, mice_format$delta, col = color_type, 
#        ylim = c(0,1), xlab = 'Freq. (Hz)', ylab = "Power Proportion",
#        main = paste0("Delta band power highlighting state", state_names[state_ind]))
#   lines(mice_format$t1, mice_format$delta, lwd = 0.1)
#   
#   plot(mice_format$t1, mice_format$theta, col = color_type, 
#        ylim = c(0,1), xlab = 'Freq. (Hz)', ylab = "Power Proportion",
#        main = paste0("Theta band power highlighting state", state_names[state_ind]))
#   lines(mice_format$t1, mice_format$theta, lwd = 0.1)
#   
#   plot(mice_format$t1, mice_format$alpha, col = color_type, 
#        ylim = c(0,1), xlab = 'Freq. (Hz)', ylab = "Power Proportion",
#        main = paste0("Alpha band power highlighting state", state_names[state_ind]))
#   lines(mice_format$t1, mice_format$alpha, lwd = 0.1)
#   
#   plot(mice_format$t1, mice_format$beta, col = color_type, 
#        ylim = c(0,1), xlab = 'Freq. (Hz)', ylab = "Power Proportion",
#        main = paste0("Beta band power highlighting state", state_names[state_ind]))
#   lines(mice_format$t1, mice_format$beta, lwd = 0.1)
#   
#   state_ind = state_ind + 1
# }
# dev.off()

pdf("stateChangeBand.pdf")
par(mfrow = c(3,1))
state_ind = 2
for(i in 1:4) {
  print(i)
  
  color_type = mice_format$state
  color_type[which(mice_format$state != 2)] = 0
  
  plot(mice_format$t1, mice_format[, i+3], col = color_type, 
       ylim = c(0,1), xlab = 'Freq. (Hz)', ylab = "Power Proportion",
       main = paste0(band_type[i], " band power highlighting state: ", state_names[2]))
  lines(mice_format$t1, mice_format[, i+3], lwd = 0.1)
  
  color_type = mice_format$state
  color_type[which(mice_format$state != 3)] = 0
  
  plot(mice_format$t1, mice_format[,i+3], col = color_type, 
       ylim = c(0,1), xlab = 'Freq. (Hz)', ylab = "Power Proportion",
       main = paste0(band_type[i], " band power highlighting state: ", state_names[3]))
  lines(mice_format$t1, mice_format[, i+3], lwd = 0.1)
  
  color_type = mice_format$state
  color_type[which(mice_format$state != 4)] = 0
  
  plot(mice_format$t1, mice_format[, i+3], col = color_type, 
       ylim = c(0,1), xlab = 'Freq. (Hz)', ylab = "Power Proportion",
       main = paste0(band_type[i], " band power highlighting state: ", state_names[4]))
  lines(mice_format$t1, mice_format[, i+3], lwd = 0.1)
  
  state_ind = state_ind + 1
}
dev.off()



