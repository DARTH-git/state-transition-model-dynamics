## Illustrative stylistic 3-state model showing the array approach ## 
rm(list = ls())  # remove any variables in R's memory 
# Load the packages
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots

# initial set up
age         <- 70  # age of starting cohort
n.t         <- 30  # time horizon, number of cycles
v.age.names <- age:(age + n.t - 1) # vector with age names
v.n <- c("H", "S", "D") # vector with the 3 health states of the model:
# Healthy (H), Sick (S), Dead (D)
n.states <- length(v.n) # number of health states 

#### Generate initial set of base-case external parameters ####
# Costs
c.H   = 1000   # cost of remaining one cycle healthy 
c.S   = 3000   # cost of remaining one cycle sick 
c.D   = 0      # cost of being dead (per cycle)
# State utilities
u.H   = 1      # utility when healthy 
u.S   = 0.60   # utility when sick 
u.D   = 0      # utility when healthy 
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
m.P["H", "S"]  <- p.HS
m.P["H", "D"]  <- p.HD
# From Sick
m.P["S", "H"]  <- p.SH
m.P["S", "S"]  <- 1 - (p.SH + p.SD)
m.P["S", "D"]  <- p.SD
# From Death
m.P["D", "H"]  <- 0
m.P["D", "S"]  <- 0
m.P["D", "D"]  <- 1

#### Initial state vector ####
v.m0 <- c(H = 1, S = 0, D = 0) # initiate the vector

## Create the Markov cohort trace matrix m.M that captures the proportion of the cohort in each state at each cycle
m.M <- matrix(0, nrow = (n.t + 1), ncol = n.states, dimnames = list(0:n.t, v.n)) # initialize cohort trace matrix 
m.M[1, ] <- v.m0   # store the initial state vector

# initiate the array 
a.A <- array(0, dim = c(n.states, n.states, n.t + 1), dimnames = list(v.n, v.n, 0:n.t)) # initialize array

diag(a.A[, , 1]) <- v.m0 # store the initial state vector in the diagonal of A

m.R.costs <- m.R.effects <- matrix(NA, nrow = n.states, ncol = n.states,  dimnames = list(v.n, v.n))

# Fill in matrix for costs
# From Healthy
m.R.costs["H", "H"]  <- c.H
m.R.costs["H", "S"]  <- c.H 
m.R.costs["H", "D"]  <- c.H + ic.D
# From Sick
m.R.costs["S", "H"]  <- c.S
m.R.costs["S", "S"]  <- c.S 
m.R.costs["S", "D"]  <- c.S + ic.D
# From Death
m.R.costs["D", "H"]  <- c.D
m.R.costs["D", "S"]  <- c.D
m.R.costs["D", "D"]  <- c.D 

# Fill in matrix for effects
# From Healthy
m.R.effects["H", "H"]  <- u.H
m.R.effects["H", "S"]  <- u.H - du.HS 
m.R.effects["H", "D"]  <- u.H 
# From Sick
m.R.effects["S", "H"]  <- u.S
m.R.effects["S", "S"]  <- u.S 
m.R.effects["S", "D"]  <- u.S 
# From Death
m.R.effects["D", "H"]  <- u.D
m.R.effects["D", "S"]  <- u.D
m.R.effects["D", "D"]  <- u.D 

#### Expected QALYs and Costs per cycle for each strategy ####
a.Y.costs <- a.Y.effects <- array(0, dim = c(n.states, n.states, n.t + 1), dimnames = list(v.n, v.n, 0:n.t))

# Initialize arrays
a.Y.costs[, , 1]   <- a.A[, , 1] * m.R.costs   
a.Y.effects[, , 1] <- a.A[, , 1] * m.R.effects 

### Run the model 
for(t in 1:n.t){  # loop through the number of cycles
  # estimate the state vector for the next cycle (t + 1)
  m.M[t + 1, ] <- m.M[t, ] %*% m.P    
  a.A[, , t + 1] <- diag(m.M[t, ]) %*% m.P  # fill array A for t + 1 
  
  # element-wise-multiplication of array A with the rewards matrices
  a.Y.costs[, , t + 1]   <- a.A[, , t + 1] * m.R.costs   
  a.Y.effects[, , t + 1] <- a.A[, , t + 1] * m.R.effects 
}

# Perform a cost-effectiveness analysis 
v.costs <- rowSums(t(colSums(a.Y.costs)))    # calculate the expected costs per cycle
v.QALYs <- rowSums(t(colSums(a.Y.effects)))  # calculate the expected QALYs per cycle
TC <- sum(v.costs)                           # calculate the total expected costs
TE <- sum(v.QALYs)                           # calculate the total expected QALYS
v.results <- c(TC, TE)                       # combine the total expected costs and QALYs
names(v.results) <- c("Costs", "Effect")     # name the vector
v.results                                    # print the results  


################################################################################
### Ratio of those that transitioned from sick to dead at each cycle to those that transitioned to dead from both healthy and sick.
v.e <- numeric(n.t)            # create the vector v.e
v.e[1] <- 0                    # initiate the vector

for(t in 1:n.t){ # calculate the ratio for all time points
  v.e[t + 1] <-  a.A["S", "D", t + 1] / (a.A["H", "D", t + 1] +  a.A["S", "D", t + 1])
}
