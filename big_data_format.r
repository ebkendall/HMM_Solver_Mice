# -----------------------------------------------------------------------------

# Investigation into each of the unique states other than the five above
# 2, 5, 10, 11, 12
# (2) has no state labels

state_names <- c("limbo", "is", "nrem", "rem", "")

for(ll in c(5,10,11,12)) {
  Dir <- Sys.glob("Data_format/Exported_ephys_and_state/*.csv")[ll]
  
  mice_data = read.csv(Dir)
  
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