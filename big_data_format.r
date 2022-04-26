# -----------------------------------------------------------------------------

# Investigation into each of the unique states other than the five above
# 2, 5, 10, 11, 12
# (2) has no state labels

state_names <- c("limbo", "is", "nrem", "rem", "")

for(ll in c(5,10,11,12)) {
 Dir <- Sys.glob("Data_format/Exported_ephys_and_state/*.csv")[ll]

 mice_data = read.csv(Dir)
 mice_data$eeg = mice_data$eeg - mean(mice_data$eeg)
 mice_data$emg = mice_data$emg - mean(mice_data$emg)
 
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

 for (i in 1:length(mice_data_sep)) {
   print(unique(mice_data_sep[[i]]$state))
   mice_data_temp = mice_data_sep[[i]]

   save(mice_data_temp, file = paste0("Data_format/Big_data_format/mice_data_",
                                      ll, "_", i, ".rda"))
 }
}

# -----------------------------------------------------------------------------

# Checking the issues with the emg and eeg data and see if we can uncover the issue
# pdf('Plots/Supplement/raw_series_both.pdf')
# par(mfrow = c(2,1))
# for (ll in 1:length(Sys.glob("Data_format/Data_raw/*.csv"))) {
#     
#   Dir <- Sys.glob("Data_format/Data_raw/*.csv")[ll]
#   print(ll)
# 
#   main_title = substring(Dir, 38)
#   print(main_title)
# 
  # mice_data = read.csv(Dir)
  # print(paste0("Mean EEG: ", mean(mice_data$ecog), " , Variance: ", sd(mice_data$ecog)))
  # print(paste0("Mean EMG: ", mean(mice_data$emg), " , Variance: ", sd(mice_data$emg)))
# 
#   mean_eeg = mean(mice_data$eeg)
#   sd_eeg = sd(mice_data$eeg)
#   print(paste0('EEG --> mean: ', round(mean_eeg, 4),', sd: ', round(sd_eeg, 4)))
# 
#   mean_emg = mean(mice_data$emg)
#   sd_emg = sd(mice_data$emg)
#   print(paste0('EMG --> mean: ', round(mean_emg, 4),', sd: ', round(sd_emg, 4)))
# 
#   plot(mice_data$eeg[1:50000], type = 'l', 
#        main = main_title, ylab = '', 
#        xlab = paste0('EEG --> mean: ', round(mean_eeg, 4),
#                      ', sd: ', round(sd_eeg, 4)))
# 
#   plot(mice_data$emg[1:50000], type = 'l', col = 'red',
#        main = main_title, ylab = '', 
#        xlab = paste0('EMG --> mean: ', round(mean_emg, 4),
#                      ', sd: ', round(sd_emg, 4)))
# 
# }
# 
# dev.off()