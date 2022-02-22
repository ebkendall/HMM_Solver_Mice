# ------------------------------------------------------------------------------
# NEW STUFF  -------------------------------------------------------------------
# ------------------------------------------------------------------------------

library(deSolve)
library(expm)

Q <- function(t,beta){

  betaMat = matrix(beta, ncol = 2, byrow = F) # determine the covariates

  q1  = exp( c(1,t) %*% betaMat[1,] )  # Transition from LIMBO to IS
  q2  = exp( c(1,t) %*% betaMat[2,] )  # Transition from LIMBO to NREM
  q3  = exp( c(1,t) %*% betaMat[3,] )  # Transition from IS    to LIMBO
  q4  = exp( c(1,t) %*% betaMat[4,] )  # Transition from IS    to NREM
  q5  = exp( c(1,t) %*% betaMat[5,] )  # Transition from IS    to REM
  q6  = exp( c(1,t) %*% betaMat[6,] )  # Transition from NREM  to LIMBO
  q7  = exp( c(1,t) %*% betaMat[7,] )  # Transition from NREM  to IS
  q8  = exp( c(1,t) %*% betaMat[8,] )  # Transition from NREM  to REM
  q9  = exp( c(1,t) %*% betaMat[9,] )  # Transition from REM   to LIMBO
  q10 = exp( c(1,t) %*% betaMat[10,] ) # Transition from REM   to IS
  q11 = exp( c(1,t) %*% betaMat[10,] ) # Transition from REM   to NREM
  q12 = exp( c(1,t) %*% betaMat[11,] ) # Transition from REM   to NREM

  # qmat = matrix(c(  0,  q1,  q2,  0,
  #                   q3,   0,  q4, q5,
  #                   q6,  q7,   0, q8,
  #                   q9, q10, q11,  0),
  #               nrow = 4, byrow = T)
  qmat = matrix(c(  0,  q1,  q2,  q3,
                    q4,   0,  q5, q6,
                    q7,  q8,   0, q9,
                    q10, q11, q12,  0),
                nrow = 4, byrow = T)
  diag(qmat) = -rowSums(qmat)

  return(qmat)
}

model_t <- function(t,p,parms) {
  qmat = Q(t, parms$b)

  pmat = matrix(c(  p[1],  p[2],  p[3],  p[4],
                    p[5],  p[6],  p[7],  p[8],
                    p[9],  p[10], p[11], p[12],
                    p[13], p[14], p[15], p[16]),
                nrow = 4, byrow = T)

  # pmat = matrix(c(  p[1],  p[2],  p[3],  0,
  #                     p[4],  p[5],  p[6],  p[7],
  #                     p[8],  p[9], p[10],  p[11],
  #                     p[12], p[13], p[14], p[15]),
  #   nrow = 4, byrow = T)

  # Vectorizing the matrix multiplication row-wise
  dP = c(t(pmat %*% qmat))

  return(list(dP))

}

init_par = pars = c(c(matrix(c(-2.26568339,  0.08766060,
                               -1.22022878, -4.44888558,
                               -1.56180104, -0.08262607,
                               -2.20978996,  0.05404948,
                               -2.41222255,  0.10833734,
                               -1.9       ,  0.07      ,
                               -1.0       ,  0.1       ,
                               -1.0       ,  0.1       ,
                               -1.0       ,  0.1       ,
                               -1.0       ,  0.1       ,
                               -1.0       ,  0.1       ), ncol=2, byrow=T)),
                    c(-3.444682, -3.850148, -4.543295,
                      -3.218876, -1.321756, -3.624341,
                      -3.624341, -1.321756, -3.218876,
                      -4.543295, -3.850148, -3.444682),
                    c( -6.52842355, -6.15970066),
                    c(1, 2, 3, 4),
                    c(1, 2, 3, 4),
                    c(1, 2, 3, 4),
                    c(1, 2, 3, 4))

# par_index = list( beta=1:20, misclass=21:32, pi_logit=33:34,
#                   l_delta = 35:38, l_theta=39:42, l_alpha=43:46, l_beta=47:50)
par_index = list( beta=1:22, misclass=23:35, pi_logit=36:37,
                  l_delta = 38:41, l_theta=42:45, l_alpha=46:49, l_beta=50:53)

parms_t <- pars[par_index$beta]
p <- c( p1=1,  p2=0, p3=0, p4=0,
           p5=0,  p6=1, p7=0, p8=0,
           p9=0, p10=0,p11=1,p12=0,
           p13=0, p14=0,p15=0,p16=1)
# p <- c( p1=1,  p2=0, p3=0,
#            p4=0,  p5=1, p6=0, p7=0,
#            p8=0, p9=0,p10=1,p11=0,
#            p12=0, p13=0,p14=0,p15=1)

times <- c(5,10)
beta <- pars[par_index$beta]

out <- deSolve::ode(p, times = times,
                             func = model_t,
                             parms = list(b=beta))

P_out <- matrix(c( out[2,"p1"],  out[2,"p2"],  out[2,"p3"],  out[2,"p4"],
               out[2,"p5"],  out[2,"p6"],  out[2,"p7"],  out[2,"p8"],
               out[2,"p9"], out[2,"p10"], out[2,"p11"], out[2,"p12"],
               out[2,"p13"], out[2,"p14"], out[2,"p15"], out[2,"p16"]),
            nrow = 4, byrow = T)

betaMat_t = matrix(parms_t, ncol = 2)

t_int <- seq(5,10, 0.001)
P = diag(4)
for(i in t_int) {

  tempQ <- Q(i, betaMat_t)
  P = P %*% expm(0.001* tempQ)

}

print("deSolve"); print(P_out)
print("expm"); print(P)


# ------------------------------------------------------------------------------
# Continuous Resp  -------------------------------------------------------------
# ------------------------------------------------------------------------------

# Q <- function(t,beta){
#   
#   betaMat = matrix(beta, ncol = 3, byrow = F) # determine the covariates
#   q1  = exp( c(1,t,1) %*% betaMat[1,] )  # Transition from IS to NREM
#   q2  = exp( c(1,t,1) %*% betaMat[2,] )  # Transition from IS to REM
#   q3  = exp( c(1,t,1) %*% betaMat[3,] )  # Transition from NREM to IS
#   q4  = exp( c(1,t,1) %*% betaMat[4,] )  # Transition from NREM to REM
#   q5  = exp( c(1,t,1) %*% betaMat[5,] )  # Transition from REM to IS
#   q6  = exp( c(1,t,1) %*% betaMat[6,] )  # Transition from REM to NREM
#   
#   qmat = matrix(c( 0,q1,q2,
#                    q3, 0,q4,
#                    q5,q6, 0),nrow=3,byrow=TRUE)
#   diag(qmat) = -rowSums(qmat)
#   
#   return(qmat)
# }
# 
# model_t <- function(t,p,parms) {
#   
#   betaMat <- matrix(parms$b, ncol = 3, byrow = F)
#   
#   q1  = exp( c(1,t,1) %*% betaMat[1,] )  # Transition from IS to NREM
#   q2  = exp( c(1,t,1) %*% betaMat[2,] )  # Transition from IS to REM
#   q3  = exp( c(1,t,1) %*% betaMat[3,] )  # Transition from NREM to IS
#   q4  = exp( c(1,t,1) %*% betaMat[4,] )  # Transition from NREM to REM
#   q5  = exp( c(1,t,1) %*% betaMat[5,] )  # Transition from REM to IS
#   q6  = exp( c(1,t,1) %*% betaMat[6,] )# Transition from REM to NREM
#   
#   dP = rep(1,9) # this is the vector with all diffqs
#   
#   dP[1] = -p[1]*(q1+q2) + p[2]*q3 + p[3]*q5
#   dP[2] = p[1]*q1 - p[2]*(q3+q4) + p[3]*q6
#   dP[3] = p[1]*q2 + p[2]*q4 - p[3]*(q5+q6)
#   dP[4] = -p[4]*(q1+q2) + p[5]*q3 + p[6]*q5
#   dP[5] = p[4]*q1 - p[5]*(q3+q4) + p[6]*q6
#   dP[6] = p[4]*q2 + p[5]*q4 - p[6]*(q5+q6)
#   dP[7] = -p[7]*(q1+q2) + p[8]*q3 + p[9]*q5
#   dP[8] = p[7]*q1 - p[8]*(q3+q4) + p[9]*q6
#   dP[9] = p[7]*q2 + p[8]*q4 - p[9]*(q5+q6)
#   
#   return(list(dP))
#   
# }
# 
# pars = init_par = trueValues = c(c(matrix(c(-2.26568339,  0.08766060, -0.49991746,
#                                      -1.22022878, -4.44888558, -0.82779213,
#                                      -1.56180104, -0.08262607,  0.73838829,
#                                      -2.20978996,  0.05404948, -1.83682627,
#                                      -2.41222255,  0.10833734,  1.63135439,
#                                      -1.9       ,  0.07      , -1.1  ), ncol=3, byrow=T)),
#                           c(  -5.73343061, -2.140066, -2.833213, -2.12144526),
#                           c( -6.52842355, -6.15970066),
#                           c(10, 20, 30),
#                           1)
# 
# par_index = list( beta=1:18, misclass=19:22, pi_logit=23:24, mu = 25:27, sigma = 28)
# 
# parms_t <- pars[par_index$beta]
# beta <- pars[par_index$beta]
# p <- c(p1=1,p2=0,p3=0,p4=0,p5=1,p6=0,p7=0,p8=0,p9=1)
# 
# times <- c(5,10)
# 
# out <- deSolve::ode(p, times = times, func = model_t, parms = list(b=beta, x_ik = 1))
# 
# P_out <- matrix(c(out[2,"p1"], out[2,"p2"], out[2,"p3"],
#               out[2,"p4"], out[2,"p5"], out[2,"p6"],
#               out[2,"p7"], out[2,"p8"], out[2,"p9"]), nrow = 3, byrow = T)
# 
# t_int <- seq(5,10, 0.001)
# P = diag(3)
# for(i in t_int) {
#   
#   tempQ <- Q(i, beta)
#   P = P %*% expm(0.001* tempQ)
#   
# }
# 
# print("deSolve"); print(P_out)
# print("expm"); print(P)
# 



# # **************************************************************************
# # This works for representing time inhomogeneous as well as time homogeneous
# # **************************************************************************
# 
# 
# Q <- function(iyears,sex,betaMat){
  
#   q1  = exp( c(1,iyears,sex) %*% betaMat[1,] )  # Transition from state 1 to state 2.
#   q2  = exp( c(1,iyears,sex) %*% betaMat[2,] )  # Transition from state 2 to state 3.
#   q3  = exp( c(1,iyears,sex) %*% betaMat[3,] )  # Transition from state 1 to death.
#   q4  = exp( c(1,iyears,sex) %*% betaMat[4,] )  # Transition from state 2 to death.
#   q5  = exp( c(1,iyears,sex) %*% betaMat[5,] )  # Transition from state 3 to death.
  
#   qmat = matrix(c( 0,q1, 0,q2,
#                    0, 0,q3,q4,
#                    0, 0, 0,q5,
#                    0, 0, 0, 0),nrow=4,byrow=TRUE) 
#   diag(qmat) = -rowSums(qmat)
  
#   return(qmat)
# }

# out_mat <- function(t, out) {
#   myMat <- matrix(c(out[t,"p1"], out[t,"p2"], out[t,"p3"], out[t,"p4"],
#                     0, out[t,"p5"], out[t,"p6"], out[t,"p7"],
#                     0,  0, out[t,"p8"], out[t,"p9"],
#                     0,  0,  0,  1), nrow = 4, byrow = T)
#   return(myMat)
# }


# #out_mat = Vectorize(out_mat, vectorize.args = "t", SIMPLIFY = T)

# # Model dependent on time ----------------------------------------------

# model_t <- function(t,p,b) {
  
#   betaMat <- matrix(b, ncol = 3, byrow = F)
  
#   q1  = exp( c(1,t,1) %*% betaMat[1,] )  # Transition from state 1 to state 2.
#   q2  = exp( c(1,t,1) %*% betaMat[2,] )  # Transition from state 2 to state 3.
#   q3  = exp( c(1,t,1) %*% betaMat[3,] )  # Transition from state 1 to death.
#   q4  = exp( c(1,t,1) %*% betaMat[4,] )  # Transition from state 2 to death.
#   q5  = exp( c(1,t,1) %*% betaMat[5,] )  # Transition from state 3 to death.
  
#   dP = rep(1,9) # this is the vector with all diffEqs
  
#   dP[1] = p[1]*(-q1-q2)
#   dP[2] = p[1]*q1 + p[2]*(-q3-q4)
#   dP[3] = p[2]*q3 - p[3]*q5
#   dP[4] = p[1]*q2 + p[2]*q4 + p[3]*q5
#   dP[5] = p[5]*(-q3-q4)
#   dP[6] = p[5]*q3 - p[6]*q5
#   dP[7] = p[5]*q4 + p[6]*q5
#   dP[8] = -p[8]*q5
#   dP[9] = p[8]*q5
  
#   return(list(dP)) 
  
# }

# trueValues <- c(c(matrix(c(-2.54, 0.11, -0.56,
#                            -2.94, -0.24,  0.15,
#                            -1.10, -0.15, -0.03,
#                            -3.92, 0.23,  0.21,
#                            -2.12, 0.08,  1.17), ncol=3, byrow=T)),
#                 c(  -4.59512, -1.15268, -2.751535, -2.090741),
#                 c( -3.178054, -4.59512))

# par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

# parms_t <- trueValues[par_index$beta]
# p <- c(p1=1,p2=0,p3=0,p4=0,p5=1,p6=0,p7=0,p8=1,p9=0)
# times <- c(5,10)

# out_time <- ode(p, times, model_t, parms = parms_t)

# P_out = out_mat(2, out_time)

# betaMat_t = matrix(parms_t, ncol = 3)

# t_int <- seq(5,10, 0.001)
# P = diag(4)
# for(i in t_int) {
  
#   tempQ <- Q(i, 1, betaMat_t)
#   P = P %*% expm(0.001* tempQ)
  
# }

# print("deSolve"); print(P_out)
# print("expm"); print(P)
