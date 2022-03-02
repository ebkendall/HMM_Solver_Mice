source("mcmc_routine.r")

# ind = as.numeric(Sys.getenv('LSB_JOBINDEX'))
args = commandArgs(TRUE)
ind = as.numeric(args[1]) # either 5 or 30 for now

set.seed(ind)

# ONLY COVARIATE IS TIME
init_par = pars = c(c(matrix(c(2.071653   , 2.85576,
                               1.212561   , -3.418799,
                               -4.198075  , 0.5761976,
                               -1.882976  , 4.516193,
                               -0.01241128, -2.303688,
                               -0.6022949 , 3.432778,
                               -0.9710429 , -2.282282,
                               1.126578   , 4.04868,
                               -0.6175124 , 0.6705703,
                               3.312069   , 0.6658912,
                               -2.792293  , 0.1919586), ncol=2, byrow=T)),
                               c(-9.266642, 0.9177677, -4.659014,
                                 -4.381079, 2.809255, -2.815062,
                                 -1.498909, -2.702002, -2.519476,
                                 -5.202432, -2.634818, -1.023371),
                               c(-2.573828,  1.353468), 
                               c(2.406984, 7.39707, 0.006235556, 0.9011201),
                               c(2.526422, 2.422841, 1.793587, 0.3687743),
                               c(-1.534735, 0.01483018, 1.649471, 5.28734),
                               c(3.082263, 2.136501, 1.184794, -0.081426944)) 

par_index = list( beta=1:22, misclass=23:34, pi_logit=35:36, 
                  l_delta = 37:40, l_theta=41:44, l_alpha=45:48, l_beta=49:52)

prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

# load("Data_format/mice_format_1.rda")
load("Data_format/mice_format.rda")

temp_data = as.matrix(mice_format); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y_1 = temp_data[,"state"]
y_2 = temp_data[,c("delta", "theta", "alpha", "beta"), drop=F]
t = temp_data[,"t1"]
steps = 20000
burnin = 5000
n_cores = 4

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("Model_out/mcmc_out_1_20_", ind, ".rda"))

# missclass_resp_fnc <- matrix(c(0.94, 0.03, 0.02, 0.01,
#                                0.03, 0.75, 0.20, 0.02,
#                                0.02, 0.20, 0.75, 0.03,
#                                0.01, 0.02, 0.03, 0.94), nrow = 4, byrow = T)

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
