source("mcmc_routine_expm.r")

# ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
ind = 3
set.seed(ind)
trialNum = 45

init_par = c(c(matrix(c(-6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0,
                        -6, 0), ncol=2, byrow = T)),
            c(-5, -5, -5, -5, -5, -5),
            c(0, 0, 0),
            log(c(0.42844057, 0.52802631, 0.46033976, 0.52947184,
              0.34666468, 0.30778351, 0.36753584, 0.30819960,
              0.17758807, 0.12382649, 0.13143568, 0.11718367,
              0.04730668, 0.04036369, 0.04068871, 0.04514489)))

par_index = list( beta=1:24, misclass = 25:30, pi_logit=31:33, l_delta = 34:37, 
                  l_theta=38:41, l_alpha=42:45, l_beta=46:49)

prior_mean = c(c(matrix(c(-8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0,
                          -8, 0), ncol=2, byrow = T)),
                    c(0, 0, 0, 0, 0, 0),
                    c(0, 0, 0),
                    c(0, 0, 0, 0, 
                      0, 0, 0, 0, 
                      0, 0, 0, 0,
                      0, 0, 0, 0))

prior_sd = c(c(matrix(c(2, 20,
                        2, 20,
                        2, 20,
                        2, 20,
                        2, 20,
                        2, 20,
                        2, 20,
                        2, 20,
                        2, 20,
                        2, 20,
                        2, 20,
                        2, 20), ncol=2, byrow = T)), 
                  c(20, 20, 20, 20, 20, 20),
                  c(20, 20, 20),
                  c(20, 20, 20, 20, 
                    20, 20, 20, 20, 
                    20, 20, 20, 20,
                    20, 20, 20, 20))


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
steps = 30000
burnin = 5000
n_cores = 2

n_post = 5000; steps2 = 30000
index_post = (steps2 - burnin - n_post + 1):(steps2 - burnin)
load(paste0('Model_out/mcmc_out_', 1, '_', trialNum - 2,'.rda'))
in_post_temp = tail(index_post, 300)
par_temp = colMeans(mcmc_out$chain[in_post_temp,])
rownames(par_temp) = NULL
par_index_old = list( beta=1:24, pi_logit=25:27, l_delta = 28:31, l_theta=32:35, 
                  l_alpha=36:39, l_beta=40:43)

init_par[par_index$beta] = par_temp[par_index_old$beta]
init_par[par_index$pi_logit] = par_temp[par_index_old$pi_logit]
init_par[par_index$l_delta] = par_temp[par_index_old$l_delta]
init_par[par_index$l_theta] = par_temp[par_index_old$l_theta]
init_par[par_index$l_alpha] = par_temp[par_index_old$l_alpha]
init_par[par_index$l_beta] = par_temp[par_index_old$l_beta]

pars = init_par

print(pars)

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores, ind, trialNum)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = 
              paste0("Model_out/mcmc_out_", ind, "_", trialNum, "_expm.rda"))
