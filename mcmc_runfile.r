source("mcmc_routine.r")

# ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
ind = 5
set.seed(ind)
trialNum = 1

# ONLY COVARIATE IS TIME (seed 6 looked like it was running the most)
init_par = c(10.7546,   0.5908,  10.3479,  -3.1937,  -3.4115,   1.8609,   5.5669,  -9.9735,
           -9.4079, -79.8430,  -6.1392,   4.5375,  -3.5097,   4.4836,   2.2844,  -0.4092,
           -7.0872,  -1.4680,  -1.8114,   8.2248,   5.6371,   3.1865,  10.7522,   1.0814,
           -6.6501, -21.0068, -18.9590, -14.3941, -18.0109, -26.0016, -21.3028, -24.4576,
          -10.1085, -18.5643,   1.1161,  -5.8191,   3.3676,   2.5901,   1.6991,   0.7753,
            2.5628,   2.7438,   1.7583,   0.8327,   1.8636,  -1.3147,  -0.8495,  -2.1251,
           -1.7250,   2.1832,  -2.9005,  -5.8171)

 
par_index = list( beta=1:22, misclass=23:34, pi_logit=35:36, 
                  l_delta = 37:40, l_theta=41:44, l_alpha=45:48, l_beta=49:52)

prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

load("Data_format/mice_format_fourMice.rda")

temp_data = as.matrix(mice_format); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y_1 = temp_data[,"state"]
y_2 = temp_data[,c("delta", "theta", "alpha", "beta"), drop=F]
t = temp_data[,"t1"]
steps = 10000
burnin = 5000
n_cores = 16

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = 
              paste0("Model_out/mcmc_out_", ind, "_", trialNum, ".rda"))


# logit_fnc <- function(x) {
#     logit_temp = matrix(0, nrow = nrow(x), ncol = ncol(x))

#     for(diag_ind in 1:nrow(x)) {
#         temp = x[diag_ind, ]

#         for(i in 1:length(temp)) {
#             if(i != diag_ind & temp[i] != 0) {
#                 logit_temp[diag_ind, i] = 
#                         log( temp[i] / (1 - sum(na.exclude(temp[-diag_ind]))))
#             }
#         }   
#     }

#     return(logit_temp)
    
# }


# Investigating how frequently the time between states changes
 # freq_df = data.frame("start" = c(1), "end" = c(1), "state" = c(1))
 # freq_df[1,] = c(mice_data$t[1], mice_data$t[1], mice_data$state[1])
 # row_num = 1
 # for (i in 2:nrow(mice_data)) {
 #   if(mice_data$state[i] != mice_data$state[i-1]) {
 #     freq_df$end[row_num] = mice_data$t[i-1]
 #     row_num = row_num + 1
 #     freq_df[row_num,] = c(mice_data$t[i], mice_data$t[i], mice_data$state[i])
 #   } else {
 #     freq_df$end[row_num] = mice_data$t[i]
 #   }
 # }

# freq_df$start = as.numeric(freq_df$start)
# freq_df$end = as.numeric(freq_df$end)
