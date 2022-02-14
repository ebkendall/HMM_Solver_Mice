file_names = c("Mice_Data/WT 08 20210309 03 Penetrating Arteriole 064 ECoG, EMG and sleep.csv",
               "Mice_Data/WT 08 20210311 06 Penetrating Arteriole 065 ECoG, EMG and sleep.csv",
               "Mice_Data/WT 08 20210311 10 Penetrating Arteriole 066 ECoG, EMG and sleep.csv",
               "Mice_Data/WT 09 20210308 03 Penetrating Arteriole 073 ECoG, EMG and sleep.csv")

# There exists 4 states at the moment and we are assuming all are possible
# To start, lets try getting rid of the "undefined" column and seeing what kind
# of transitions we get
for (i in 1:4) {
  mice_data = read.csv(file_names[i])
  
  # for(j in 1:nrow(mice_data)) {
  #   if(mice_data$state[j] != "<undefined>") {
  #     print(paste0("Time: ", mice_data$t[j], ", State: ", mice_data$state[j]))
  #   }
  # }
  
  transition_mat = matrix(0, nrow = 4, ncol = 4)
  colnames(transition_mat) = c("<undefined>", "IS", "NREM", "REM")
  rownames(transition_mat) = c("<undefined>", "IS", "NREM", "REM")

  for (i in 1:(nrow(mice_data) - 1)) {
    s1 = mice_data$state[i]
    s2 = mice_data$state[i+1]
    row_num = which(rownames(transition_mat) == s1)
    col_num = which(colnames(transition_mat) == s2)
    transition_mat[row_num, col_num] = transition_mat[row_num, col_num] + 1
  }

  mice_data = mice_data[mice_data$state != "<undefined>", ]

  transition_mat2 = matrix(0, nrow = 3, ncol = 3)
  colnames(transition_mat2) = c("IS", "NREM", "REM")
  rownames(transition_mat2) = c("IS", "NREM", "REM")

  for (i in 1:(nrow(mice_data) - 1)) {
    s1 = mice_data$state[i]
    s2 = mice_data$state[i+1]
    row_num = which(rownames(transition_mat2) == s1)
    col_num = which(colnames(transition_mat2) == s2)
    transition_mat2[row_num, col_num] = transition_mat2[row_num, col_num] + 1
  }

  print("All States")
  print(transition_mat)
  print("No undefined")
  print(transition_mat2)
}

