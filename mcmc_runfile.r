source("mcmc_routine.r")

# ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
ind = 5
set.seed(ind)
trialNum = 1

# ONLY COVARIATE IS TIME 
# init_par = pars =  c(3.39430821,  10.12405236,  -7.93689935,   4.07974000,   2.08356686,
#                     -9.52969326,   3.27907507, -15.91484866,  -9.28285114,  -4.56532190,
#                     1.49610183,  -4.91409187, -12.84450383,  -1.76026395,   0.03354449,
#                     0.07772913,  -0.41705965,   0.61153827,   0.94644509,   0.04946497,
#                     0.34727929,   0.71425303,  -2.27113849,  -2.85133668,  -1.45095693,
#                     -14.00300984, -17.11449746,  -4.83359798, -12.57249206,  -9.05233237,
#                     -15.16377749, -12.04371566, -11.79666308, -11.91959835,   6.33681539,
#                     5.53072276,
#              log(8 * c(0.45301142, 0.29606903, 0.15807648, 0.09284307,
#              0.3159233,  0.3298007,  0.2507649,  0.1035110,
#              0.43535440, 0.30149348, 0.17704457, 0.08610755,
#              0.41606167, 0.36792410, 0.15595449, 0.06005974)))

init_par = pars = c(10.71203032,   3.52090244,  11.10473555,   0.36116910,   2.57988950,
 -5.08623890,   4.23518359,  -0.65745901,   3.65054620, -81.04296032,
 -9.97705498,   4.47413475,  -0.08061642,   4.37178174,   0.51065118,
 -0.16417233,  -5.50421397,   0.02489742,  -0.51344842,  -0.13201315,
 -2.45150251,   1.28375292,  13.87753704,  -0.03992211,  11.98177040,
 -5.53781652, -18.84627889, -28.98020087, -10.43218983, -22.07223472,
-24.04738660, -11.78338951,  -3.43488772, -10.88070673,  -0.99765845,
 -5.15287242,   3.54751261,   3.22850625,   2.27816601,   1.29591907,
  2.84414351,   3.31423511,   2.07661581,   1.19061879,   3.51267867,
  2.40335638,   1.49102911,   0.63561830,   2.26214127,   2.27664348,
  1.98225601,   0.85473975)

par_index = list( beta=1:22, misclass=23:34, pi_logit=35:36, 
                  l_delta = 37:40, l_theta=41:44, l_alpha=45:48, l_beta=49:52)

prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

# load('Data_format/mice_format_sub_total.rda')
load('Data_format/mice_format_sub_total_15.rda')
mice_format = mice_format[!(mice_format$ptnum %in% c(18,19,20,21,22,23)), ]

temp_data = as.matrix(mice_format); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y_1 = temp_data[,"state"]
y_2 = temp_data[,c("delta", "theta", "alpha", "beta"), drop=F]
t = temp_data[,"t1"]
steps = 10000
burnin = 5000
n_cores = 5

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = 
              paste0("Model_out/mcmc_out_", ind, "_", trialNum, ".rda"))
