source("mcmc_routine.r")

ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(ind)
trialNum = 25

init_par = pars =  c(10.0050,   4.1675,  12.1816,  -1.9492,   3.6943,  -3.1640,   3.9072,  -2.9780,
                      2.2770, -77.7631, -14.1756,   3.5643,   0.4594,   3.0818,   0.8064,  -0.6416,
                     -3.1180,   0.5931,  -0.1604,   0.3956,  -1.1921,   4.6106,  13.2611,   1.6711,
                      7.6776,  -5.9033, -17.9331, -32.6448,  -7.4676, -24.5494, -25.3553, -10.4891,
                     -4.3325, -12.3260,   0.1338,  -4.7456,   3.8262,   3.3905,   2.4872,   1.5384,
                      2.7699,   3.0584,   2.2187,   1.0079,   3.4377,   2.6117,   1.6008,   0.7604,
                      2.5420,   2.1901,   1.2825,   0.6821)

par_index = list( beta=1:22, misclass=23:34, pi_logit=35:36, 
                  l_delta = 37:40, l_theta=41:44, l_alpha=45:48, l_beta=49:52)

prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

load('Data_format/mice_format_sub_total_split.rda')
mice_format = mice_format[!(mice_format$ptnum %in% 34:59), ]

sub_ind = c(101,  77,   3,  97,  70,  84, 113,  82,  30,  88,   2,  24,   1,  61,  27,  74,  99, 104,  32,
             73, 110,  15,  22,  90,  98,  75,  80,  11,  16,  14,  21,  78, 102,  91, 103,  96,  13,  20,
             12,  29,   5,  67,  68,  85,   8,  10,  60, 106,  92,  87)
 
mice_format = mice_format[(mice_format$ptnum %in% sub_ind), ]

print(table(mice_format$state))

temp_data = as.matrix(mice_format); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y_1 = temp_data[,"state"]
y_2 = temp_data[,c("delta", "theta", "alpha", "beta"), drop=F]
t = temp_data[,"t1"]
steps = 25000
burnin = 5000
n_cores = 20

n_post = 5000; steps2 = 25000
index_post = (steps2 - burnin - n_post + 1):(steps2 - burnin)
load(paste0('Model_out/mcmc_out_', 4, '_', trialNum - 1,'.rda'))
in_post_temp = tail(index_post, 300)
par_temp = colMeans(mcmc_out$chain[in_post_temp,])
rownames(par_temp) = NULL
init_par = par_temp
pars = par_temp

print(pars)

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores, ind, trialNum)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = 
              paste0("Model_out/mcmc_out_", ind, "_", trialNum, ".rda"))
