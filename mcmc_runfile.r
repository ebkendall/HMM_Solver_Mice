source("mcmc_routine.r")

#ind = as.numeric(Sys.getenv('LSB_JOBINDEX'))
ind=10
set.seed(ind)

# ONLY COVARIATE IS TIME (seed 6 looked like it was running the most)
init_par = pars = c(c(matrix(c(-3.5795145 , 1.8837027,
                               -1.6320443 , -0.01087329,
                               -2.54012846 , -0.41003774,
                               1.7573955 , 0.69776741,
                               2.2151124 , 0.36636279,
                               2.77588470 , -0.88555272,
                               -1.0840123 , 1.258807628,
                               -1.2000797 , 0.8709151,
                               1.2516498 , 0.008028053,
                               -0.8830492 , 2.55136088,
                               -0.97033314 , 2.4799403), ncol=2, byrow=T)),
                               c(-2.369395, -3.318072, -4.821319,
                                 -4.116913, -2.012214, -5.524204,
                                 -3.264101, -4.335017, -3.851207,
                                 -0.8467074, -4.250335, -0.4214324),
                               c(-0.52303533, -0.72652087), 
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
steps = 50000
burnin = 5000
n_cores = 16

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = 
              paste0("Model_out/mcmc_out_", ind, "_", steps / 1000, "_30.rda"))


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

# init_par = pars = c(c(matrix(c(2.071653   , 2.85576,
#                                1.212561   , -3.418799,
#                                -4.198075  , 0.5761976,
#                                -1.882976  , 4.516193,
#                                -0.01241128, -2.303688,
#                                -0.6022949 , 3.432778,
#                                -0.9710429 , -2.282282,
#                                1.126578   , 4.04868,
#                                -0.6175124 , 0.6705703,
#                                3.312069   , 0.6658912,
#                                -2.792293  , 0.1919586), ncol=2, byrow=T)),
#                                c(-9.266642, 0.9177677, -4.659014,
#                                  -4.381079, 2.809255, -2.815062,
#                                  -1.498909, -2.702002, -2.519476,
#                                  -5.202432, -2.634818, -1.023371),
#                                c(-2.573828,  1.353468), 
#                                c(2.406984, 7.39707, 0.006235556, 0.9011201),
#                                c(2.526422, 2.422841, 1.793587, 0.3687743),
#                                c(-1.534735, 0.01483018, 1.649471, 5.28734),
#                                c(3.082263, 2.136501, 1.184794, -0.081426944)) 
