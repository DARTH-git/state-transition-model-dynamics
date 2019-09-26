## Illustrative stylistic 3-state model showing the tranditional appraoch ## 

# Appendix of 'Krijkamp EM, Alarid-Escudero F, Enns EA, Hunink MGM, Pechlivanoglou P, Jalal HJ.' **A multidimensional array representation of state-transition model dynamics.** (revised, September 2019)

# Load the packages
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots

# initial set up
age         <- 70  # age of starting cohort
n_t         <- 30  # time horizon, number of cycles
v_age_names <- age:(age + n_t - 1) # vector with age names
v_n <- c("H", "Stemp", "S", "Dtemp", "D") # vector with the 3 health states of the model:
# Healthy (H), Sick (S), Dead (D) and two temporary health states one for Sick for the first time (Stemp) and one for dying (Dtemp)
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
m_P["H", "H"]      <- 1 - (p_HS + p_HD)
m_P["H", "Stemp"]  <- p_HS
m_P["H", "S"]      <- 0
m_P["H", "Dtemp"]  <- p_HD
m_P["H", "D"]      <- 0

# From Sick temporary (first cycle being sick)
m_P["Stemp", "H"]      <- p_SH
m_P["Stemp", "Stemp"]  <- 0
m_P["Stemp", "S"]      <- 1 - (p_SH + p_SD)
m_P["Stemp", "Dtemp"]  <- p_SD
m_P["Stemp", "D"]      <- 0

# From Sick
m_P["S", "H"]      <- p_SH
m_P["S", "Stemp"]  <- 0
m_P["S", "S"]      <- 1 - (p_SH + p_SD)
m_P["S", "Dtemp"]  <- p_SD
m_P["S", "D"]      <- 0

# From Death temporary
m_P["Dtemp", "H"]     <- 0
m_P["Dtemp", "Stemp"] <- 0
m_P["Dtemp", "S"]     <- 0
m_P["Dtemp", "Dtemp"] <- 0
m_P["Dtemp", "D"]     <- 1

# From Death
m_P["D", "H"]      <- 0
m_P["D", "Stemp"]  <- 0
m_P["D", "S"]      <- 0
m_P["D", "Dtemp"]  <- 0
m_P["D", "D"]      <- 1

#### Cohort trace matrix ####
## Initial state vector
v_m0 <- c(H = 1, Stemp = 0, S = 0, D = 0, Dtemp = 0) # all the cohort starts in the Healthy state

## Create the Markov cohort trace matrix m_M that captures the proportion of 
## the cohort in each state at each cycle
m_M <- matrix(0, 
              nrow = (n_t + 1), 
              ncol = n_states, 
              dimnames = list(0:n_t, v_n)) # initialize cohort trace matrix 
m_M[1, ] <- v_m0 # store the initial state vector in the first row of the cohort trace


#### State and tranisition rewards ####
## Create a vector to store rewards
v_R_costs   <- c(c_H, c_S, c_S, c_D + ic_D, c_D)
v_R_effects <- c(u_H, u_S - du_HS, u_S, u_D, u_D)
names(v_R_costs) <- names(v_R_effects) <- v_n

#### Expected QALYs and Costs per cycle for each strategy ####
## Create matrix to store expected outcomes
m_Y_costs <- m_Y_effects <- matrix(0, 
                                   nrow = n_states, 
                                   ncol = n_states,  
                                   dimnames = list(v_n, v_n))



#### Run the cSTM ####
for(t in 1:n_t){  # loop through the number of cycles
  # estimate the state vector for the next cycle (t + 1)
  m_M[t + 1, ] <- m_M[t, ] %*% m_P    
}

#### Aggregate outcomes ####
v_costs <- m_M %*% v_R_costs                # calculate the expected costs per cycle
v_QALYs <- m_M %*% v_R_effects              # calculate the expected QALYs per cycle

# apply upfront costs
v_costs_upfront <- v_effect_upfront <- c(0, 0, 0, 0, 0)
v_costs[1] <- m_M[1, ] %*% v_costs_upfront 
v_QALYs[1] <- m_M[1, ] %*% v_effect_upfront

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

