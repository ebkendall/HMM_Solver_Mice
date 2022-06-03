source("mcmc_routine.r")

ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(ind)
trialNum = 25

init_par = pars =  c(c(matrix(c(3.31784266, -0.03402560,
                                2.83361407, -0.08131234,
                                2.83361407, -0.08131234, # new row
                                3.84319607, -0.02059774,
                                0.65920888, 0.05320999,
                                -0.10363139, 0.54404285,
                                3.93103199, -0.05186608, 
                                1.64904538, 0.24675486,
                                -11.93696963, -8.67149268,
                                0.96464474, 0.25204027,
                                -44.11083783, 3.88994039,
                                -9.95628008, 2.54984249), ncol=2, byrow = T)),
                    c(-0.50031735, -30.89456518, -18.57665356, -17.15469869,
                      -6.95465942, -8.53797330), 
                    c(1.25816427, -1.15710610),
                    c(3.53763678, 2.89439678, 1.95569138, 1.04666927, 
                      2.61727157, 2.65737401, 1.81185671, 0.82234054, 
                      3.90673024, 2.54431988, 1.53650704, 0.73327927,
                      1.90583103, 1.77204698, 0.72612163, 0.04022586))

par_index = list( beta=1:24, misclass=25:30, pi_logit=31:32, 
                  l_delta = 33:36, l_theta=37:40, l_alpha=41:44, l_beta=45:48)

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

# n_post = 5000; steps2 = 25000
# index_post = (steps2 - burnin - n_post + 1):(steps2 - burnin)
# load(paste0('Model_out/mcmc_out_', 4, '_', trialNum - 1,'.rda'))
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
