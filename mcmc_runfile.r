source("mcmc_routine.r")

# ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
ind = 5
set.seed(ind)
trialNum = 1

# ONLY COVARIATE IS TIME 
init_par = c(-1.30120910,   2.44636060,  -4.73585806,   2.84463848,   2.97213573,  -6.56580193, 
              1.72444643,  -5.26770296,  -6.58936483,  -3.84855576,   1.35129114,   0.05878274, 
             -3.34609834,  -1.62550665,  -0.37374388,  -0.37447001,  -0.67107971,   0.88727980, 
              0.93825006,  -0.42982539,   0.21440851,   0.86510679,  -2.86740732,  -2.51391481, 
             -2.86123630,  -9.06509746,  -2.61245196,  -2.13697775, -10.44760980,  -2.77309567, 
            -14.73438731, -10.16502261,  -8.84188090,  -9.84476672,   5.02655881,   4.26384594,
             rep(1, 16))

            
par_index = list( beta=1:22, misclass=23:34, pi_logit=35:36, 
                  l_delta = 37:40, l_theta=41:44, l_alpha=45:48, l_beta=49:52)

prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

load("Data_format/mice_format_total_new.rda")

temp_data = as.matrix(mice_format); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y_1 = temp_data[,"state"]
y_2 = temp_data[,c("delta", "theta", "alpha", "beta"), drop=F]
t = temp_data[,"t1"]
steps = 20000
burnin = 5000
n_cores = 16

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = 
              paste0("Model_out/mcmc_out_", ind, "_", trialNum, ".rda"))
