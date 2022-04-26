library(tidyverse)
library(gridExtra)
library(grid)

# load("../Data_format/mice_format_ecog.rda")
load('Data_format/mice_format_total_new_5sec.rda')
mice_format = mice_format[mice_format$ptnum == 2, ]
df = unique(diff(mice_format$t1))[1]

mice_format$state = as.factor(mice_format$state)
# levels(mice_format$state) = c("IS", "NREM", "REM", "undefined")

x_s = (max(mice_format$t1) - min(mice_format$t1)) / 8
x_m = min(mice_format$t1)

# eeg_plot = vector(mode = "list", length = 2)
# eeg_plot[[1]] = 

ggplot(mice_format, aes(t1, state)) + 
  geom_rect(aes(NULL, NULL, 
                xmin=t1, xmax=t1+df, 
                ymin=0, ymax=1, 
                fill=state)) + 
  scale_fill_manual(values = alpha(c("blue", "red", "green", "yellow"), .15)) +
  geom_line(aes(y = delta), color = "blue", linetype = "solid") +
  geom_line(aes(y = theta), color = "red", linetype = "solid") +
  geom_line(aes(y = alpha), color = "purple", linetype = "solid") +
  geom_line(aes(y = beta), color = "black", linetype = "solid") +
  guides(color = guide_legend(order=1)) + 
  xlab("Time (center & scale)") +
  ylab("Proportion of power") +
  ggtitle("Relationship between state and band powers (EEG)") +
  theme(text = element_text(size = 20)) +
  annotation_custom(textGrob(' -- Delta -- ', gp = gpar(col = "blue", fontface = 'bold')), 
                    xmin = x_m + 2*x_s, xmax = x_m + 2*x_s, ymin = 1.02, ymax = 1.02) +
  annotation_custom(textGrob(' -- Theta -- ', gp = gpar(col = "red", fontface = 'bold')), 
                    xmin = x_m + 3*x_s, xmax = x_m + 3*x_s, ymin = 1.02, ymax = 1.02) +
  annotation_custom(textGrob(' -- Alpha -- ', gp = gpar(col = "purple", fontface = 'bold')), 
                    xmin = x_m + 4*x_s, xmax = x_m + 4*x_s, ymin = 1.02, ymax = 1.02) +
  annotation_custom(textGrob(' -- Beta -- ', gp = gpar(col = 'black', fontface = 'bold')), 
                    xmin = x_m + 5*x_s, xmax = x_m + 5*x_s, ymin = 1.02, ymax = 1.02)

# load("../Data_format/mice_format_emg.rda")
# library(tidyverse)
# 
# mice_format$state = as.factor(mice_format$state)
# levels(mice_format$state) = c("IS", "NREM", "REM", "undefined")
# 
# eeg_plot[[2]] = ggplot(mice_format, aes(t1, state)) + 
#   geom_rect(aes(NULL, NULL, 
#                 xmin=t1, xmax=t1+5, 
#                 ymin=0, ymax=1, 
#                 fill=state)) + 
#   scale_fill_manual(values = alpha(c("blue", "red", "green", "yellow"), .15)) +
#   geom_line(aes(y = delta), color = "blue", linetype = "solid") +
#   geom_line(aes(y = theta), color = "red", linetype = "solid") +
#   geom_line(aes(y = alpha), color = "purple", linetype = "solid") +
#   geom_line(aes(y = beta), color = "black", linetype = "solid") +
#   guides(color = guide_legend(order=1)) + 
#   xlab("Time") +
#   ylab("Proportion of power") +
#   ggtitle("Relationship between state and band powers (EMG)") +
#   theme(text = element_text(size = 20)) +
#   annotation_custom(textGrob(' -- Delta -- ', gp = gpar(col = "blue", fontface = 'bold')), 
#                     xmin = 250, xmax = 250, ymin = 1.02, ymax = 1.02) +
#   annotation_custom(textGrob(' -- Theta -- ', gp = gpar(col = "red", fontface = 'bold')), 
#                     xmin = 500, xmax = 500, ymin = 1.02, ymax = 1.02) +
#   annotation_custom(textGrob(' -- Alpha -- ', gp = gpar(col = "purple", fontface = 'bold')), 
#                     xmin = 750, xmax = 750, ymin = 1.02, ymax = 1.02) +
#   annotation_custom(textGrob(' -- Beta -- ', gp = gpar(col = 'black', fontface = 'bold')), 
#                     xmin = 1000, xmax = 1000, ymin = 1.02, ymax = 1.02)
# 
# ggsave(eeg_plot[[1]], filename = "../Plots/ecog.pdf", height = 8.5,
#   width = 11)
# ggsave(eeg_plot[[2]], filename = "../Plots/emg.pdf", height = 8.5,
#   width = 11)

# How often do we observe state changes ----------------------------------------

# mice_data = read.csv("~/Dropbox/Shared_HMM_ICU/mouse_data/WT 08 20210309 03 Penetrating Arteriole 064 ECoG, EMG and sleep.csv")

# state_change = data.frame("start" = mice_data$t[1], 
#                           "end" = mice_data$t[1],
#                           "state" = mice_data$state[1])
# ind = 1

# for(i in 2:nrow(mice_data)) {
#   if(mice_data$state[i] != mice_data$state[i-1]) {
#     state_change[ind, ] = c(state_change$start[ind], 
#                             mice_data$t[i], 
#                             state_change$state[ind])
#     ind = ind + 1
#     state_change[ind, ] = c(mice_data$t[i], mice_data$t[i], mice_data$state[i])
#   } else {
#     state_change$end[ind] = mice_data$t[i]
#   }
# }
