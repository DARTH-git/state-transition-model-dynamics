# This code is used to prove the cumulative incidence calculations based on the 
# array structure by running a simple 3-state model:
# (Healthy (H), Sick (S), Dead(D))
# with an additional health state for those that have always been healthy (H) and 
# one for those that have been sick at least once before (H2)
# The cumulative incidence is the:
# sum of the cohort that transitioned from H -> S / total individuals in H at t=0


## Simple 3-state model showing the array approach ## 
rm(list = ls())  # remove any variables in R's memory 
# Load the packages
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots

# initial set up
age         <- 70  # age of starting cohort
n.t         <- 30  # time horizon, number of cycles
v.age.names <- age:(age + n.t - 1) # vector with age names
v.n <- c("H", "H2", "S", "D") # vector with the 3 health states of the model:
# Healthy (H), Sick (S), Dead (D)
n.states <- length(v.n) # number of health states 

# Transition probabilities (per cycle)
p.HS  = 0.30   # probability to become sick when healthy
p.HD  = 0.05   # probability to die when healthy
p.SH  = 0.15   # probability to become healthy when sick
p.SD  = 0.20   # probability to die when healthy
# Transition rewards
du.HS  = 0.10  # one-time utility decrement when becoming sick
ic.D   = 4000   # one-time cost of dying

#### Transition probability matrix ####
# matrix m.P at the first cycle
m.P <- matrix(NA, nrow = n.states, ncol = n.states, dimnames = list(v.n, v.n))

# Fill in matrix
# From Healthy
m.P["H", "H"]  <- 1 - (p.HS + p.HD)
m.P["H", "H2"] <- 0
m.P["H", "S"]  <- p.HS
m.P["H", "D"]  <- p.HD
# From Healthy2
m.P["H2", "H"]  <- 0
m.P["H2", "H2"] <- 1 - (p.HS + p.HD)
m.P["H2", "S"]  <- p.HS
m.P["H2", "D"]  <- p.HD
# From Sick
m.P["S", "H"]   <- 0
m.P["S", "H2"]  <- p.SH
m.P["S", "S"]   <- 1 - (p.SH + p.SD)
m.P["S", "D"]   <- p.SD
# From Death
m.P["D", "H"]  <- 0
m.P["D", "H2"]  <- 0
m.P["D", "S"]  <- 0
m.P["D", "D"]  <- 1

#### Initial state vector ####
v.m0 <- c(H = 1, H2 = 0,  S = 0, D = 0) # initiate the vector

## Create the Markov cohort trace matrix m.M that capturs the proportion of the cohort in each state at each cycle
m.M <- matrix(0, nrow = (n.t + 1), ncol = n.states, dimnames = list(0:n.t, v.n)) # initialize cohort trace matrix 
m.M[1, ] <- v.m0   # store the initial state vector

# initiate the array 
a.A <- array(0, dim = c(n.states, n.states, n.t + 1), dimnames = list(v.n, v.n, 0:n.t)) # initialize array

diag(a.A[, , 1]) <- v.m0 # store the initial state vector in the diagonal of A

m.R.costs <- m.R.effects <- matrix(NA, nrow = n.states, ncol = n.states,  dimnames = list(v.n, v.n))


### Run the model 
for(t in 1:n.t){  # loop through the number of cycles
  # estimate the state vector for the next cycle (t + 1)
  m.M[t + 1, ] <- m.M[t, ] %*% m.P    
  a.A[, , t + 1] <- diag(m.M[t, ]) %*% m.P  # fill array A for t + 1 
  
}


### Plot cohort trace
ggplot(melt(m.M), aes(x = Var1, y = value, color = Var2)) +
  geom_line(size = 1.3) +
  scale_color_discrete(l = 50, name = "Health state", h = c(45, 365)) +
  xlab("Cycle") +
  ylab("Proportion of the cohort") +
  ggtitle("Cohort trace of the simple 3-state model")+
  theme_bw(base_size = 16) +
  scale_x_continuous(name = "Cycles", limits = c(0, n.t), breaks = seq(0, n.t, 10)) +
  theme()


### Cumulative incidence 

# the cumulative incidence are all of those that made the transition form healthy (H) to sick over time.
# this is stored in the a.A["H", "S", ] data
# those at risk of becomming a new case are those that stay in the H state
# this information if valuable and therefore stored in a .csv file to compare with the results of the 3-state model.

m.M_4states <- m.M # store the values in a different object name 
a.A_4states <- a.A # store the values in a different object name 
  
# Save the objects
save(m.M_4states, file = "output/Cohort_trace_4states.RData") # save the object
save(a.A_4states, file = "output/Array_4states.RData") # save the object 





