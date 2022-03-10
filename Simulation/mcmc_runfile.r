source("mcmc_routine.r")

ind = as.numeric(Sys.getenv('LSB_JOBINDEX'))

set.seed(ind)

# ONLY COVARIATE IS TIME (seed 6 looked like it was running the most)
init_par = pars = c(c(matrix(c(-2.26568339,  0.08766060,
                      -1.22022878, -4.44888558,
                      -1.56180104, -0.08262607,
                      -2.20978996,  0.05404948,
                      -2.41222255,  0.10833734,
                      -1.9       ,  0.07      ,
                      -1.0       ,  0.1       ,
                      -1.0       ,  0.1       ,
                      -1.0       ,  0.1       ,
                      -1.0       ,  0.1       ,
                      -1.0       ,  0.1       ), ncol=2, byrow=T)),
           c(-3.444682, -3.850148, -4.543295,
             -3.218876, -1.321756, -3.624341,
             -3.624341, -1.321756, -3.218876,
             -4.543295, -3.850148, -3.444682),
           c( -6.52842355, -6.15970066),
           c(1, 1, 1, 1),
           c(1.25, 1, 0.5, 0),
           c(1.5, 0.8, 0.5, 0),
           c(1, 1, 0.5, 0))

par_index = list( beta=1:22, misclass=23:34, pi_logit=35:36, 
                  l_delta = 37:40, l_theta=41:44, l_alpha=45:48, l_beta=49:52)

prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

load(paste0("Data_simulation/miceData", ind, ".rda"))

temp_data = as.matrix(miceData); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y_1 = temp_data[,"state"]
y_2 = temp_data[,c("delta", "theta", "alpha", "beta"), drop=F]
t = temp_data[,"time"]
steps = 20000
burnin = 5000
n_cores = 16

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = 
              paste0("Model_out/mcmc_out_", ind, "_", steps / 1000, ".rda"))


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