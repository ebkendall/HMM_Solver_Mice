source("mcmc_routine.r")

# ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
ind=10
set.seed(ind)

# ONLY COVARIATE IS TIME (seed 6 looked like it was running the most)
init_par = c(1.85061099,   1.27008193,  -9.13595201,   1.60538168,  -0.14610977,
       -8.96542692,   1.77124043, -11.27758079,  -9.96261824, -87.77478076,
        0.59369285,  -4.37142690,  -1.93922920,  -0.19847169,  -0.14707156,
        0.12525060,  -0.34632909,  -0.23345521,  -0.28437221,  -3.83891689,
        6.54804750,  -0.08665573,  -2.86256107,  -2.91963863,  -2.98420832,
      -12.77693040, -25.37430065, -12.74349187, -12.88851272, -23.16303582,
      -12.85471183, -12.43304076, -12.44444766, -12.46997750,  -0.52300000,
       -0.72650000,
     c(0, 0, 0, 0),
     c(0, 0, 0, 0),
     c(0, 0, 0, 0),
     c(0, 0, 0, 0))



par_index = list( beta=1:22, misclass=23:34, pi_logit=35:36, 
                  l_delta = 37:40, l_theta=41:44, l_alpha=45:48, l_beta=49:52)

prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

load("Data_format/mice_format_total.rda")

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
              paste0("Model_out/mcmc_out_", ind, "_", steps / 1000, "_5.rda"))


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
