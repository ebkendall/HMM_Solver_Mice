library(signal, quietly = T)
library(eegkit, quietly = T)

# ------------------------------------------------------------------------------
# Simulated data to get a better understanding of what this band pass filter is doing
# ------------------------------------------------------------------------------

t = seq(0, 0.5, length.out = 500)
f = 10* sin(2*pi*t*50) +  20*sin(2 * pi * t * 150)  + 30* sin(2*pi*t*250)
plot(t, f, type = "l")

omega = seq(0, 1/(t[2] - t[1]), length.out =length(t))

fft_sin = fft(f)

eeg_filt_100 = eegfilter(f, Fs = 1/(t[2] - t[1]), lower = 1, upper = 100, method = "butter")
eeg_filt_200 = eegfilter(f, Fs = 1/(t[2] - t[1]), lower = 100, upper = 200, method = "butter")
eeg_filt_300 = eegfilter(f, Fs = 1/(t[2] - t[1]), lower = 200, upper = 300, method = "butter")

plot(t, eeg_filt_300, type = "l", col = "blue")
lines(t, eeg_filt_100, col = "red")
lines(t, eeg_filt_200, col = "green")

fft_eeg_100 = fft(eeg_filt_100)
fft_eeg_200 = fft(eeg_filt_200)
fft_eeg_300 = fft(eeg_filt_300)

plot(omega, abs(fft_eeg_100) / length(fft_eeg_100), type = "h")
plot(omega, abs(fft_eeg_200) / length(fft_eeg_200), type = "h")
plot(omega, abs(fft_eeg_300) / length(fft_eeg_300), type = "h")

# Illustrates if we can reform the original plot
plot(t, f, type = "l")
lines(t, eeg_filt_100 + eeg_filt_200 + eeg_filt_300, col = "green")

# ------------------------------------------------------------------------------
# Investigating the difference between our method and the built in function
# ------------------------------------------------------------------------------
omega_color = omega
omega_color[which(omega < 100)] = 1
omega_color[which(omega >= 100 & omega < 200)] = 2
omega_color[which(omega >= 200 & omega < 300)] = 3
omega_color[which(omega >= 300)] = 4
plot(omega[1:250], abs(fft_sin)[1:250] / length(fft_sin), col = omega_color,
     type = "h", xlab = "HZ")

plot(omega[1:250], abs(fft_eeg_100)[1:250] / length(fft_eeg_100), 
     type = "h", ylim = c(0,15),xlim = c(0,500), col = 1, xlab= "", ylab = "")
par(new=TRUE)
plot(omega[1:250], abs(fft_eeg_200)[1:250] / length(fft_eeg_200), 
     type = "h", ylim = c(0,15),xlim = c(0,500), col = 2, xlab= "", ylab = "")
par(new=TRUE)
plot(omega[1:250], abs(fft_eeg_300)[1:250] / length(fft_eeg_300), 
     type = "h", ylim = c(0,15),xlim = c(0,500), col = 3, xlab= "", ylab = "")

# Lets see how these two different plots relate to one another
# First we will look at the red coloring
fft_red_orig = fft_sin[which(omega >= 100 & omega < 200)]
fft_red_fnc  = fft(eeg_filt_200)

pow_red_orig = sum(abs(fft_red_orig)) / (2 * pi * length(fft_red_orig))
pow_red_fnc  = sum(abs(fft_red_fnc)) / (2 * pi * length(fft_red_fnc))

print(paste0("Power for 100 - 200 Hz using our method: ", pow_red_orig))
print(paste0("Power for 100 - 200 Hz using eegfilter method: ", pow_red_fnc))


# ------------------------------------------------------------------------------
# Doing the above illustrationg with the real data
# ------------------------------------------------------------------------------

Dir <- "~/Dropbox/Shared_HMM_ICU/mouse_data/WT 08 20210309 03 Penetrating Arteriole 064 ECoG, EMG and sleep.csv"

mice_data = read.csv(Dir)

index_seq = seq(1, nrow(mice_data), 2560)

i=1

test = mice_data[index_seq[i]:index_seq[i+1], ] # grabs data every 5 seconds

plot(test$t, test$ecog, type = "l")

# f_delta = eegfilter(test$ecog, Fs = 512, 
#                     lower = 0.001, upper = 4, method = "butter")
# f_theta = eegfilter(test$ecog, Fs = 512, 
#                     lower = 4, upper = 7, method = "butter")
# f_alpha = eegfilter(test$ecog, Fs = 512, 
#                     lower = 7, upper = 12, method = "butter")
# f_beta  = eegfilter(test$ecog, Fs = 512, 
#                     lower = 12, upper = 30, method = "butter")
# 
# plot(test$t, test$ecog, type = "l", lty = 1, lwd = 2, ylim = c(-1, 1))
# lines(test$t, f_delta, col = "red", lty = 2, lwd = 2)
# lines(test$t, f_theta, col = "blue", lty = 2, lwd = 2)
# lines(test$t, f_alpha, col = "green", lty = 2, lwd = 2)
# lines(test$t, f_beta, col = "purple", lty = 2, lwd = 2)
# 
# 
# fft_delta = fft(f_delta)
# fft_theta = fft(f_theta)
# fft_alpha = fft(f_alpha)
# fft_beta  = fft(f_beta)

# omega = seq(0, 1/(test$t[2] - test$t[1]), length.out =length(test$t))
# omega_color = omega
# omega_color[which(omega < 4)] = 1
# omega_color[which(omega >= 4 & omega < 7)] = 2
# omega_color[which(omega >= 7 & omega < 12)] = 3
# omega_color[which(omega >= 12 & omega < 30)] = 4
# omega_color[which(omega >= 30)] = 5

fft_new = fft(test$ecog)

pdf("fft_band_pass.pdf")
par(mfrow = c(2,1))

plot(omega[1:500], abs(fft_new)[1:500] / length(fft_new), 
     type = "h", ylim = c(0,0.041), col = omega_color, xlab= "", ylab = "",
     main = "FFT of the original EEG signal", xlim = c(0,40))
legend(80, 0.04, legend=c("< 4 Hz", "4 - 7 Hz", "7 - 12 Hz", "12 - 30 Hz", "> 30 Hz"),
       col=c(1,2,3,4,5),
       title="Hertz Filter", fill=1:5, text.font=4)


plot(omega[1:500], abs(fft_delta)[1:500] / length(fft_delta), 
     type = "h", ylim = c(0,0.041),xlim = c(0,40), col = 1, xlab= "", ylab = "")
par(new=TRUE)
plot(omega[1:500], abs(fft_theta)[1:500] / length(fft_theta), 
     type = "h", ylim = c(0,0.041),xlim = c(0,40), col = 2, xlab= "", ylab = "")
par(new=TRUE)
plot(omega[1:500], abs(fft_alpha)[1:500] / length(fft_alpha), 
     type = "h", ylim = c(0,0.041),xlim = c(0,40), col = 3, xlab= "", ylab = "")
par(new=TRUE)
plot(omega[1:500], abs(fft_beta)[1:500] / length(fft_beta),
     type = "h", ylim = c(0,0.041),xlim = c(0,40), col = 4, xlab= "", ylab = "")
legend(80, 0.04, legend=c("< 4 Hz", "4 - 7 Hz", "7 - 12 Hz", "12 - 30 Hz"),
       col=c(1,2,3,4),
       title="Hertz Filter", fill=1:4, text.font=4)

# Illustrates if we can reform the original plot
plot(test$t, test$ecog, type = "l", main = "Recovering the original time series")
lines(test$t, f_delta + f_theta + f_alpha + f_beta, col = "green")

dev.off()
