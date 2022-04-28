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