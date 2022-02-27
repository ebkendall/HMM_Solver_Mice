# ------------------------------------------------------------------------------
# Simulated data to get a better understanding of what this band pass filter is doing
# ------------------------------------------------------------------------------

t = seq(0, 0.5, length.out = 500)
f = sin(1200 * 2 * pi * t) +  0.5 * sin(90 * 2 * pi * t) # + 30* sin(2*pi*t*250)
plot(t, f, type = "l")

omega = seq(0, 1/(t[2] - t[1]), length.out =length(t))

fft_sin = fft(f)

# ------------------------------------------------------------------------------
# Investigating the difference between our method and the built in function
# ------------------------------------------------------------------------------
plot(omega, abs(fft_sin) / length(fft_sin), type = "h", xlab = "HZ", ylab = "")



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

K = nrow(test)
k = 1:K
w_k = 0.5 - 0.5 * cos(2 * pi * k / (K-1))
hann_rw = test$ecog * w_k

plot(test$t, hann_rw, type = "l")

omega = seq(0, 1/(test$t[2] - test$t[1]), length.out =length(test$t))
omega_color = omega
omega_color[which(omega < 4)] = 1
omega_color[which(omega >= 4 & omega < 7)] = 2
omega_color[which(omega >= 7 & omega < 12)] = 3
omega_color[which(omega >= 12 & omega < 30)] = 4
omega_color[which(omega >= 30)] = 5

fft_new = fft(test$ecog)
fft_han = fft(hann_rw)

pdf("fft_compare_hann.pdf")
par(mfrow = c(2,1))

plot(omega[1:500], abs(fft_new)[1:500] / length(fft_new), 
     type = "h", ylim = c(0,0.041), col = omega_color, xlab= "", ylab = "",
     main = "FFT of the original EEG signal", xlim = c(0,40))
legend(30, 0.04, legend=c("< 4 Hz", "4 - 7 Hz", "7 - 12 Hz", "12 - 30 Hz", "> 30 Hz"),
       col=c(1,2,3,4,5),
       title="Hertz Filter", fill=1:5, text.font=4)
plot(omega[1:500], abs(fft_han)[1:500] / length(fft_han), 
     type = "h", ylim = c(0,0.041), col = omega_color, xlab= "", ylab = "",
     main = "FFT of the original EEG signal", xlim = c(0,40))
legend(30, 0.04, legend=c("< 4 Hz", "4 - 7 Hz", "7 - 12 Hz", "12 - 30 Hz", "> 30 Hz"),
       col=c(1,2,3,4,5),
       title="Hertz Filter", fill=1:5, text.font=4)

dev.off()
