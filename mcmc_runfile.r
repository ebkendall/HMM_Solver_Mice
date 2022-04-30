source("mcmc_routine.r")

ind = as.numeric(Sys.getenv('LSB_JOBINDEX'))
set.seed(ind)
trialNum = 12

init_par = pars = c(10.9526,   3.5101,  11.3124,  -0.5893,   2.9171,  -3.2827,   4.2641,  -0.5597,
                     4.0112, -77.8469, -14.5411,   4.2361,  -0.0658,   4.1675,  -0.3550,  -0.0404,
                    -4.4737,   0.0192,  -0.9096,  -0.0582,  -1.4610,   4.0299,  11.9184,   1.5325,
                     6.5380,  -6.5788, -16.2629, -31.6463,  -8.0170, -25.3148, -24.1587, -13.0393,
                    -5.2219, -12.4362,   0.2193,  -4.8980,   3.5589,   3.2329,   2.2567,   1.2796,
                     2.8662,   3.3651,   2.1222,   1.2287,   3.4878,   2.3617,   1.4663,   0.6263,
                     2.3741,   2.3018,   2.1155,   0.9213)

par_index = list( beta=1:22, misclass=23:34, pi_logit=35:36, 
                  l_delta = 37:40, l_theta=41:44, l_alpha=45:48, l_beta=49:52)

prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

load('Data_format/mice_format_sub_total_15.rda')
mice_format = mice_format[!(mice_format$ptnum %in% c(18,19,20,21,22,23)), ]

temp_data = as.matrix(mice_format); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y_1 = temp_data[,"state"]
y_2 = temp_data[,c("delta", "theta", "alpha", "beta"), drop=F]
t = temp_data[,"t1"]
steps = 10000
burnin = 5000
n_cores = 15

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = 
              paste0("Model_out/mcmc_out_", ind, "_", trialNum, ".rda"))
