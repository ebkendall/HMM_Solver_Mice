# Downsampling the data

for (ll in 3:length(Sys.glob("Data_format/Exported_ephys_and_state/*.csv"))) {
  
  Dir <- Sys.glob("Data_format/Exported_ephys_and_state/*.csv")[ll]
  print(paste0(ll, ", ", Dir))
  
  mice_data_big = read.csv(Dir)
  
  mice_ind = seq(1, nrow(mice_data_big), 9)
  
  mice_data = mice_data_big[mice_ind, ]
  
  save(mice_data, file = paste0('Data_format/Exported_ephys_downsample/mice_data_',
                                ll, '.rda'))

}

# -----------------------------------------------------------------------------

# Investigation into each of the unique states other than the five above

state_names <- c("limbo", "is", "nrem", "rem", "")

for(ll in c(1,3:12)) {

   load(paste0('Data_format/Exported_ephys_downsample/mice_data_', ll, '.rda'))

   change_ind = which(mice_data$state != dplyr::lag(mice_data$state))
  
   mice_data_sep = vector(mode = 'list', length = 1)
   ind = 1; start_ind = 1
  
   for(i in 1:length(change_ind)) {
     if(!(mice_data$state[change_ind[i]] %in% state_names)) {
       print(paste0(mice_data$state[change_ind[i] - 1], " --> ", mice_data$state[change_ind[i]]))
       mice_data_sep[[ind]] = mice_data[start_ind:(change_ind[i] - 1), ]
       start_ind = change_ind[i+1]
       ind = ind + 1
     }
   }
    
   if(!is.null(mice_data_sep[[1]])) {
     for (i in 1:length(mice_data_sep)) {
       print(unique(mice_data_sep[[i]]$state))
       mice_data_small = mice_data_sep[[i]]
       
       save(mice_data_small, file = paste0("Data_format/Exported_ephys_downsample/mice_data_",
                                          ll, "_", i, ".rda"))
     }
   }
   
}

# -----------------------------------------------------------------------------

# editing the series so that we have more subjects with frequency of 556

for(ll in 1:length(Sys.glob("Data_format/Exported_ephys_downsample/*.rda"))) {
  print(ll)
  
  Dir <- Sys.glob("Data_format/Exported_ephys_downsample/*.rda")[ll]
  
  load(Dir)
  
  # For the ones that had to be broken up
  small_check = substring(Dir, 39)
  if(nchar(small_check) > 16) mice_data = mice_data_small
  
  mice_data_temp = mice_data
  
  length_CONST = 614381
  
  if(nrow(mice_data_temp) <= length_CONST) {
    mice_data = mice_data_temp
    save(mice_data, file = paste0("Data_format/Exported_ephys_downsample_2/mice_data_",
                                  ll, "_2.rda"))
  } else {
    ind_breakdown = c(seq(1, nrow(mice_data_temp), length_CONST), nrow(mice_data_temp) + 1)
    
    for(ii in 2:length(ind_breakdown)) {
      mice_ind = (ind_breakdown[ii - 1]):(ind_breakdown[ii] - 1)
      mice_data = mice_data_temp[mice_ind, ]
      save(mice_data, file = paste0("Data_format/Exported_ephys_downsample_2/mice_data_",
                                    ll, "_", ii, ".rda"))
    }
  }
}



# -----------------------------------------------------------------------------

# Checking the issues with the emg and eeg data and see if we can uncover the issue
pdf('Plots/Supplement/raw_series_both_sub.pdf')
par(mfrow = c(2,1))
for (ll in 1:length(Sys.glob("Data_format/Exported_ephys_downsample/*.rda"))) {

  Dir <- Sys.glob("Data_format/Exported_ephys_downsample/*.rda")[ll]
  print(ll)

  main_title = Dir
  print(main_title)

  load(Dir)
  print(paste0("Mean EEG: ", mean(mice_data$eeg), " , Variance: ", sd(mice_data$eeg)))
  print(paste0("Mean EMG: ", mean(mice_data$emg), " , Variance: ", sd(mice_data$emg)))

  mean_eeg = mean(mice_data$eeg)
  sd_eeg = sd(mice_data$eeg)
  print(paste0('EEG --> mean: ', round(mean_eeg, 4),', sd: ', round(sd_eeg, 4)))

  mean_emg = mean(mice_data$emg)
  sd_emg = sd(mice_data$emg)
  print(paste0('EMG --> mean: ', round(mean_emg, 4),', sd: ', round(sd_emg, 4)))

  plot(mice_data$eeg[1:50000], type = 'l',
       main = main_title, ylab = '',
       xlab = paste0('EEG --> mean: ', round(mean_eeg, 4),
                     ', sd: ', round(sd_eeg, 4)))

  plot(mice_data$emg[1:50000], type = 'l', col = 'red',
       main = main_title, ylab = '',
       xlab = paste0('EMG --> mean: ', round(mean_emg, 4),
                     ', sd: ', round(sd_emg, 4)))

}

dev.off()

# -----------------------------------------------------------------------------

# New histograms
load('Data_format/mice_format_sub_total_split.rda')
# load('Data_format/mice_format_sub_total_15.rda')
# Looking into the 5 second distribution
temp = mice_format[mice_format$delta < 0.2, ]
table(temp$ptnum)
table(mice_format$ptnum)
# t_full = as.data.frame(table(mice_format$ptnum))
# t_sub = as.data.frame(table(temp$ptnum))
# t_sub$Freq / t_full$Freq
# which((t_sub$Freq / t_full$Freq) < 0.9) #c(1,3,4,5,6)
# 
mice_format = mice_format[!(mice_format$ptnum %in% 34:59), ]
# mice_format = mice_format[mice_format$ptnum == 1, ]
# 
pdf('Plots/Supplement/mice_format_sanityCheck.pdf')
par(mfrow = c(2,2))
hist(mice_format$delta[mice_format$state == 1],
     main = "Delta, state 2")
hist(mice_format$delta[mice_format$state == 2],
     main = "Delta, state 3")
hist(mice_format$delta[mice_format$state == 3],
     main = "Delta, state 4")
hist(mice_format$delta[mice_format$state == 99],
     main = "Delta, state 99")
hist(mice_format$theta[mice_format$state == 1],
     main = "Theta, state 2")
hist(mice_format$theta[mice_format$state == 2],
     main = "Theta, state 3")
hist(mice_format$theta[mice_format$state == 3],
     main = "Theta, state 4")
hist(mice_format$theta[mice_format$state == 99],
     main = "Theta, state 99")
hist(mice_format$alpha[mice_format$state == 1],
     main = "alpha, state 2")
hist(mice_format$alpha[mice_format$state == 2],
     main = "alpha, state 3")
hist(mice_format$alpha[mice_format$state == 3],
     main = "alpha, state 4")
hist(mice_format$alpha[mice_format$state == 99],
     main = "alpha, state 99")
hist(mice_format$beta[mice_format$state == 1],
     main = "beta, state 2")
hist(mice_format$beta[mice_format$state == 2],
     main = "beta, state 3")
hist(mice_format$beta[mice_format$state == 3],
     main = "beta, state 4")
hist(mice_format$beta[mice_format$state == 99],
     main = "beta, state 99")

dev.off()
