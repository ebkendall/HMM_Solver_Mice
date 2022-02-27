source("mcmc_routine.r")

# ind = as.numeric(Sys.getenv('LSB_JOBINDEX'))
ind = 1

set.seed(ind)

# ONLY COVARIATE IS TIME
init_par = pars = c(c(matrix(c(-1 , 0.05,
                               -1 , 0.05,
                               -1 , 0.05,
                               -1 , 0.05,
                               -1 , 0.05,
                               -1 , 0.05,
                               -1 , 0.05,
                               -1 , 0.05,
                               -1 , 0.05,
                               -1 , 0.05,
                               -1 , 0.05), ncol=2, byrow=T)),
                               c(-4, -4, -4,
                                 -4, -4, -4,
                                 -4, -4, -4,
                                 -4, -4, -4),
                               c(0.074049718,  0.075399602), 
                               c(2.1, 2.2, 2.3, 2.4),
                               c(2.1, 2.2, 2.3, 2.4),
                               c(2.1, 2.2, 2.3, 2.4),
                               c(2.1, 2.2, 2.3, 2.4)) 

par_index = list( beta=1:22, misclass=23:34, pi_logit=35:36, 
                  l_delta = 37:40, l_theta=41:44, l_alpha=45:48, l_beta=49:52)

prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

# load("Data_format/mice_format_1.rda")
load("Data_format/mice_format_center.rda")

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

save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, ".rda"))

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
