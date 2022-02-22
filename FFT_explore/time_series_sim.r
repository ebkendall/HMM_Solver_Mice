# ----------------------------------------
# ------- Initialization & Notes ---------
# ----------------------------------------

# State 1: a_1 = 0, a_2 = 1     (pure beta rhythm)
# State 2: a_1 = 0.5, a_2 = 0.5 (mixture rhythm)
# State 3: a_1 = 0, a_2 = 1     (pure alpha rhythm)

# *** The initial state is 1 ***

A = matrix(c(0  ,   1,
             0.5, 0.5,
             0  ,   1), nrow = 3, byrow = T)

rownames(A) = c("State_1", "State_2", "State_3")


f = c(10, 25)

P = matrix(c(0.8, 0.2, 0.0,
             0.0, 0.8, 0.2,
             0.2, 0.0, 0.8), nrow = 3, byrow = T)

t = seq(0, 40, by = 0.01)

# (1) Generate true state sequence (first a Markov process (not hidden))
true_state = c(1)

for(i in 2:length(t)) {
    true_state[i] = sample(1:3, size=1, prob=P[true_state[i-1],]) 
}

# (2) Generate the response function
y = rep(0, length(true_state))

for(i in 1:length(t)) {
  
  s = true_state[i]
  
  y[i] = A[s,1] * sin(2*pi*f[1]*t[i]) + A[s,2] * sin(2*pi*f[2]*t[i]) + rnorm(1, 0, 0.3)
  
}

AR_coeff <-c()

t_start = t_end = 0

while(t_end < 40 - 0.4) {
  print(t_start)
  
  t_end = t_start + 0.5
  
  t_sub = t[which(t <= t_end & t >= t_start)]
  y_sub = y[which(t <= t_end & t >= t_start)]
  
  temp = arima(y_sub, c(8,0,0))
  
  AR_coeff = cbind(AR_coeff, temp$coef[1:8])
  
  t_start = t_end - 0.1
  
}


