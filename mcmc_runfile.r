source("mcmc_routine_new.r")

ind = as.numeric(Sys.getenv('LSB_JOBINDEX'))

set.seed(ind)

# ONLY COVARIATE IS TIME
init_par = pars = c(c(matrix(c(-2.26568339,  0.08766060,
                               -1.22022878, -4.44888558,
                               -1.56180104, -0.08262607,
                               -2.20978996,  0.05404948,
                               -2.41222255,  0.10833734,
                               -1.9       ,  0.07      ,
                               -1.0       ,  0.1       ,
                               -1.0       ,  0.1       ,
                               -1.0       ,  0.1       ,
                               -1.0       ,  0.1       ), ncol=2, byrow=T)),
                               c(-3.444682, -3.850148, -4.543295,
                                 -3.218876, -1.321756, -3.624341,
                                 -3.624341, -1.321756, -3.218876,
                                 -4.543295, -3.850148, -3.444682),
                               c( -6.52842355, -6.15970066),
                               c(10, 20, 30, 40),
                               c(10, 20, 30, 40),
                               c(10, 20, 30, 40),
                               c(10, 20, 30, 40))

par_index = list( beta=1:20, misclass=21:32, pi_logit=33:34, 
                  l_delta = 35:38, l_theta=39:42, l_alpha=43:46, l_beta=47:50)

prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

load(paste0("Data_format/mice_format_", ind, ".rda"))

temp_data = as.matrix(miceData); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y_1 = temp_data[,"state"]
y_2 = temp_data[,c("delta", "theta", "alpha", "beta"), drop=F]
x = temp_data[, c("years", "sex"), drop=F]
t = temp_data[,"years"]
steps = 10000
burnin = 5000
n_cores = 16

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, x, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("Model_out/mcmc_out_", ind, ".rda"))

# missclass_resp_fnc <- matrix(c(0.94, 0.03, 0.02, 0.01,
#                                0.03, 0.75, 0.20, 0.02,
#                                0.02, 0.20, 0.75, 0.03,
#                                0.01, 0.02, 0.03, 0.94), nrow = 4, byrow = T)
#
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