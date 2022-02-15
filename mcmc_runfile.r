source("mcmc_routine_new.r")

ind = as.numeric(Sys.getenv('LSB_JOBINDEX'))

set.seed(ind)

init_par = trueValues = c(c(matrix(c(-2.26568339,  0.08766060, -0.49991746,
                                     -1.22022878, -4.44888558, -0.82779213,
                                     -1.56180104, -0.08262607,  0.73838829,
                                     -2.20978996,  0.05404948, -1.83682627,
                                     -2.41222255,  0.10833734,  1.63135439,
                                     -1.9       ,  0.07      , -1.1  ), ncol=3, byrow=T)),
                                     c(  -5.73343061, -0.78623894, -2.52747176, -2.12144526),
                                     c( -6.52842355, -6.15970066),
                                     c(10, 20, 30),
                                     1)

par_index = list( beta=1:18, misclass=19:22, pi_logit=23:24, mu = 25:27, sigma = 28)

betaMat <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)

prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

load(paste0("DataOut/Continuous/miceData", ind, ".rda"))

temp_data = as.matrix(miceData); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y_1 = temp_data[,"state"]
y_2 = temp_data[,"cont_resp"]
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
