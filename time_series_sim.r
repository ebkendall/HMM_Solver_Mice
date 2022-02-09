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


# (2) Generate the response function
