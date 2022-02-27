load("../Data_format/mice_format_1.rda")
library(tidyverse)

mice_format$state = as.factor(mice_format$state)

ggplot(mice_format, aes(t1, state)) + 
  geom_rect(aes(NULL, NULL, 
                xmin=t1, xmax=t1+5, 
                ymin=0, ymax=1, 
                fill=state)) + 
  scale_fill_manual(values = alpha(c("blue", "red", "green", "yellow"), .15)) +
  geom_line(aes(y = delta), color = "blue", linetype = "solid") +
  geom_line(aes(y = theta), color = "red", linetype = "solid") +
  geom_line(aes(y = alpha), color = "purple", linetype = "solid") +
  geom_line(aes(y = beta), color = "black", linetype = "solid") +
  guides(color = guide_legend(order=1)) + 
  xlab("Time") +
  ylab("Proportion of power") +
  ggtitle("Relationship between state and band powers (Hann Window)") +
  theme(text = element_text(size = 20))


# How often do we observe state changes ----------------------------------------

mice_data = read.csv("~/Dropbox/Shared_HMM_ICU/mouse_data/WT 08 20210309 03 Penetrating Arteriole 064 ECoG, EMG and sleep.csv")

state_change = data.frame("start" = mice_data$t[1], 
                          "end" = mice_data$t[1],
                          "state" = mice_data$state[1])
ind = 1

for(i in 2:nrow(mice_data)) {
  if(mice_data$state[i] != mice_data$state[i-1]) {
    state_change[ind, ] = c(state_change$start[ind], 
                            mice_data$t[i], 
                            state_change$state[ind])
    ind = ind + 1
    state_change[ind, ] = c(mice_data$t[i], mice_data$t[i], mice_data$state[i])
  } else {
    state_change$end[ind] = mice_data$t[i]
  }
}
