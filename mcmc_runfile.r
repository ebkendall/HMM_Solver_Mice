source("mcmc_routine.r")

ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(ind)
trialNum = 32

# init_par = pars =  c(c(matrix(c(3.31784266, -0.03402560,
#                                 2.83361407, -0.08131234,
#                                 2.83361407, -0.08131234, # new row
#                                 3.84319607, -0.02059774,
#                                 0.65920888, 0.05320999,
#                                -0.10363139, 0.54404285), ncol=2, byrow = T)),
#                     # c(-4.5, -4.5, -4.5, -4.5, -4.5, -4.5), 
#                     c(-1.15710610),
#                     c(3.53763678, 2.89439678, 3.90673024, 
#                       2.61727157, 2.65737401, 2.54431988, 
#                       1.77204698, 1.90583103, 1.53650704,
#                       0.72612163, 0.72612163, 0.72612163))

# par_index = list( beta=1:12, misclass=13:18, pi_logit=19, 
#                   l_delta = 20:22, l_theta=23:25, l_alpha=26:28, l_beta=29:31)


init_par = pars =  c(c(matrix(c(-1, 0,
                                -1, 0,
                                -1, 0,
                                -1, 0,
                                -1, 0,
                                -1, 0), ncol=2, byrow = T)),
                    c(0, 0),
                    log(c(0.42548916, 0.54504416, 0.34652136,
                      0.35085130, 0.30013012, 0.38799344,
                      0.16622615, 0.11803219, 0.19766568,
                      0.05743339, 0.03679353, 0.06781952)))
# Start at the sample means for 

par_index = list( beta=1:12, pi_logit=13:14, l_delta = 15:17, l_theta=18:20, 
                  l_alpha=21:23, l_beta=24:26)



prior_mean = c(c(matrix(c(0, 0,
                          0, 0,
                          0, 0, # new row
                          0, 0,
                          0, 0,
                          0, 0), ncol=2, byrow = T)),
                    # c(-4.5, -4.5, -4.5, -4.5, -4.5, -4.5), 
                    c(0, 0),
                    c(0, 0, 0, 
                      0, 0, 0, 
                      0, 0, 0,
                      0, 0, 0))

prior_sd = c(c(matrix(c(20, 20,
                        20, 20,
                        20, 20, # new row
                        20, 20,
                        20, 20,
                        20, 20), ncol=2, byrow = T)),
                  # c(7/6, 7/6, 7/6, 7/6, 7/6, 7/6), 
                  c(20, 20),
                  c(20, 20, 20, 
                    20, 20, 20, 
                    20, 20, 20,
                    20, 20, 20))


prior_par = data.frame( prior_mean= prior_mean,
                        prior_sd= prior_sd)

load('Data_format/mice_format_sub_total_split.rda')
mice_format = mice_format[!(mice_format$ptnum %in% 34:59), ]

sub_ind = c(101,  77,   3,  97,  70,  84, 113,  82,  30,  88,   2,  24,   1,  61,  27,  74,  99, 104,  32,
             73, 110,  15,  22,  90,  98,  75,  80,  11,  16,  14,  21,  78, 102,  91, 103,  96,  13,  20,
             12,  29,   5,  67,  68,  85,   8,  10,  60, 106,  92,  87)
 
mice_format = mice_format[(mice_format$ptnum %in% sub_ind), ]

temp_data = as.matrix(mice_format); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y_1 = temp_data[,"state"]
y_2 = temp_data[,c("delta", "theta", "alpha", "beta"), drop=F]
t = temp_data[,"t1"]
steps = 25000
burnin = 5000
n_cores = 2

# n_post = 5000; steps2 = 25000
# index_post = (steps2 - burnin - n_post + 1):(steps2 - burnin)
# load(paste0('Model_out/mcmc_out_', 7, '_', trialNum - 1,'.rda'))
# in_post_temp = tail(index_post, 300)
# par_temp = colMeans(mcmc_out$chain[in_post_temp,])
# rownames(par_temp) = NULL
# init_par = par_temp
# pars = par_temp

print(pars)

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores, ind, trialNum)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = 
              paste0("Model_out/mcmc_out_", ind, "_", trialNum, ".rda"))
