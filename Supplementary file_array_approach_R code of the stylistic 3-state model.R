## Illustrative stylistic 3-state model showing the array approach ## 

# Appendix of 'Krijkamp EM, Alarid-Escudero F, Enns EA, Hunink MGM, Pechlivanoglou P, Jalal HJ.' **A multidimensional array representation of state-transition model dynamics.** (revised, September 2019)

# Load the packages
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots

# initial set up
age         <- 70  # age of starting cohort
n_t         <- 30  # number of cycles
v_age_names <- age:(age + n_t - 1) # vector with age names
v_n <- c("H", "S", "D") # vector with the 3 health states of the model:
# Healthy (H), Sick (S), Dead (D)
n_states <- length(v_n) # number of health states 

#### Generate initial set of base-case external parameters ####
# Costs
c_H   <- 1000   # cost of remaining one cycle healthy 
c_S   <- 3000   # cost of remaining one cycle sick 
c_D   <- 0      # cost of being dead (per cycle)
# State utilities
u_H   <- 1      # utility when healthy 
u_S   <- 0.60   # utility when sick 
u_D   <- 0      # utility when healthy 
# Transition probabilities (per cycle)
p_HS  <- 0.30   # probability to become sick when healthy
p_HD  <- 0.05   # probability to die when healthy
p_SH  <- 0.15   # probability to become healthy when sick
p_SD  <- 0.20   # probability to die when sick
# Transition rewards
du_HS <- 0.10  # one-time utility decrement when becoming sick
ic_D  <- 4000   # one-time cost of dying

#### Transition probability matrix ####
# matrix m_P at the first cycle
m_P <- matrix(NA, 
              nrow = n_states, 
              ncol = n_states, 
              dimnames = list(v_n, v_n))

# Fill in matrix
# From Healthy
m_P["H", "H"]  <- 1 - (p_HS + p_HD)
m_P["H", "S"]  <- p_HS
m_P["H", "D"]  <- p_HD
# From Sick
m_P["S", "H"]  <- p_SH
m_P["S", "S"]  <- 1 - (p_SH + p_SD)
m_P["S", "D"]  <- p_SD
# From Death
m_P["D", "H"]  <- 0
m_P["D", "S"]  <- 0
m_P["D", "D"]  <- 1

#### Cohort trace matrix ####
## Initial state vector
v_m0 <- c(H = 1, S = 0, D = 0) # all the cohort starts in the Healthy state

## Create the Markov cohort trace matrix m_M that captures the proportion of 
## the cohort in each state at each cycle
m_M <- matrix(0, 
              nrow = (n_t + 1), 
              ncol = n_states, 
              dimnames = list(0:n_t, v_n)) # initialize cohort trace matrix 
m_M[1, ] <- v_m0 # store the initial state vector in the first row of the cohort trace

#### Multidimensional array ####
## Create the multidimensional array a_A that captures the proportion of the 
## cohort that tranistiones between health states at each cycle
a_A <- array(0, 
             dim = c(n_states, n_states, n_t + 1), 
             dimnames = list(v_n, v_n, 0:n_t)) # initialize multidimensional array

diag(a_A[, , 1]) <- v_m0 # store the initial state vector in the diagonal of the first slice of A

#### State and tranisition rewards ####
## Create matrices to store rewards
m_R_costs <- m_R_effects <- matrix(NA, 
                                   nrow = n_states, 
                                   ncol = n_states,  
                                   dimnames = list(v_n, v_n))

# Fill in matrix for costs
# To Healthy
m_R_costs["H", "H"]  <- c_H
m_R_costs["S", "H"]  <- c_H 
m_R_costs["D", "H"]  <- c_H 
# To Sick
m_R_costs["H", "S"]  <- c_S
m_R_costs["S", "S"]  <- c_S 
m_R_costs["D", "S"]  <- c_S 
# To Death
m_R_costs["H", "D"]  <- c_D + ic_D
m_R_costs["S", "D"]  <- c_D + ic_D
m_R_costs["D", "D"]  <- c_D 

# Fill in matrix for effects
# To Healthy
m_R_effects["H", "H"]  <- u_H
m_R_effects["S", "H"]  <- u_H 
m_R_effects["D", "H"]  <- u_H 
# To Sick
m_R_effects["H", "S"]  <- u_S - du_HS
m_R_effects["S", "S"]  <- u_S 
m_R_effects["D", "S"]  <- u_S 
# To Death
m_R_effects["H", "D"]  <- u_D
m_R_effects["S", "D"]  <- u_D
m_R_effects["D", "D"]  <- u_D 

#### Expected QALYs and Costs per cycle for each strategy ####
## Create multidimensional arrays to store expected outcomes
a_Y_costs <- a_Y_effects <- array(0, 
                                  dim = c(n_states, n_states, n_t + 1), 
                                  dimnames = list(v_n, v_n, 0:n_t))

# Initialize arrays
a_Y_costs[, , 1]   <- a_A[, , 1] * m_R_costs   
a_Y_effects[, , 1] <- a_A[, , 1] * m_R_effects 

#### Run the cSTM ####
for(t in 1:n_t){  # loop through the number of cycles
  # estimate the state vector for the next cycle (t + 1)
  m_M[t + 1, ] <- m_M[t, ] %*% m_P    
  a_A[, , t + 1] <- diag(m_M[t, ]) %*% m_P  # estimate the transition dynamics at t + 1
  
  # element-wise-multiplication of array A with the rewards matrices
  a_Y_costs[, , t + 1]   <- a_A[, , t + 1] * m_R_costs   
  a_Y_effects[, , t + 1] <- a_A[, , t + 1] * m_R_effects 
}

#### Aggregate outcomes ####
v_costs <- rowSums(t(colSums(a_Y_costs)))    # calculate the expected costs per cycle
v_QALYs <- rowSums(t(colSums(a_Y_effects)))  # calculate the expected QALYs per cycle
TC <- sum(v_costs)                           # calculate the total expected costs
TE <- sum(v_QALYs)                           # calculate the total expected QALYS
v_results <- c(TC, TE)                       # combine the total expected costs and QALYs
names(v_results) <- c("Costs", "Effect")     # name the vector
v_results                                    # print the results  


################################################################################
### Ratio of those that transitioned from sick to dead at each cycle to those that transitioned to dead from both healthy and sick.
v_e <- numeric(n_t + 1)   # create the vector v_e
v_e[1] <- 0               # initiate the vector

### calculate the ratio across all cycles starting in cycle 2
v_e[-1] <-  a_A["S", "D", -1] / (a_A["H", "D", -1] +  a_A["S", "D", -1])

