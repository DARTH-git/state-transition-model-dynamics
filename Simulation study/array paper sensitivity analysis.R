# Sensitivity analysis on the cohort trace approach vs. multidimensional array approach.

if (!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)

# Function to randomly sample probabilities for the transition probabilities
rand.sum <- function (n) {
  x <- sort(runif(n - 1)) # sample (n-1) random numbers from the uniform distribution and sort them
  c(x, 1) - c(0, x) # calculate the length of intervals 
}

# traditional cohort trace approach
run_trace_approach <- function (n_states, n_t) {
  # n_states: number of states
  # n_t:      number of cycles
  
  set.seed(1) 
  # total number of states (assume no absorbing states) is 
  # n_states * (n_states - 1) temporary states + the n_states original states  
  # that is, every state has (n_states - 1) temporary states 
  n_states    <- n_states * (n_states - 1) + n_states 
  m_P         <- t(replicate(n_states, rand.sum(n_states)))  # transition probability matrix 
  m_M         <- matrix(0, nrow = n_t + 1, ncol = n_states)  # initiate the Markov trace
  v_R_costs   <- matrix(rnorm(n_states, mean = 500, sd = 150), nrow = n_states, ncol = 1)  # vector of costs
  v_R_effects <- matrix(runif(n_states), nrow = n_states, ncol = 1)  # vector of QALYs
  
  #### initialize the trace 
  v_m0     <- rep(0, n_states) # everybody starts at the first state
  v_m0[1]  <- 1
  m_M[1, ] <- v_m0
  
  #### run the model
  for(t in 1:n_t){
    m_M[t + 1, ] <- m_M[t, ] %*% m_P
  }
  
  #### Aggregate outcomes 
  v_costs <- m_M %*% v_R_costs              # calculate the expected costs per cycle
  v_QALYs <- m_M %*% v_R_effects            # calculate the expected QALYs per cycle
  
  TC <- sum(v_costs)                        # calculate the total expected costs
  TE <- sum(v_QALYs)                        # calculate the total expected QALYs
  v_results <- c(TC, TE)                    # combine the total expected costs and QALYs
  names(v_results) <- c("Costs", "Effect")  # name the vector
  
  #### memory
  # store all relevant objects and compute the memory it takes
  stored <- list(m_P, m_M, v_R_costs, v_R_effects)
  mem    <- object.size(stored)
  
  output <- list(results = v_results, memory = mem)
  return(output)
}

# test function
run_trace_approach(n_states = 3, n_t = 60)


# Dynamics-array approach
run_array_approach <- function(n_states, n_t){
  # n_states: number of states
  # n_t:      number of cycles
  
  set.seed(1)
  m_P <- t(replicate(n_states, rand.sum(n_states)))     # transition probability matrix
  m_M <- matrix(0, nrow = (n_t + 1), ncol = n_states)   # initiate the Markov trace 
  a_A <- array(0, dim = c(n_states, n_states, n_t + 1)) # initiate the transition dynamics array
  
  m_R_costs   <- matrix(rnorm(n_states ^ 2, mean = 500, sd = 150), nrow = n_states, ncol = n_states)  # time-homogenous effect matrix
  m_R_effects <- matrix(runif(n_states ^ 2), nrow = n_states, ncol = n_states)                        # time-homogenous effect matrix
  a_Y_costs   <- a_Y_effects <- array(0, dim = c(n_states, n_states, n_t + 1))                        # initiate the rewards arrays
  
  #### initialize the trace and the arrays
  # everybody starts at the first state
  v_m0     <- rep(0, n_states)
  v_m0[1]  <- 1
  m_M[1, ] <- v_m0
  
  # store the initial state vector in the diagonal of the first slice of A
  diag(a_A[, , 1]) <- v_m0   
  
  # Initialize arrays
  a_Y_costs[, , 1]   <- a_A[, , 1] * m_R_costs   
  a_Y_effects[, , 1] <- a_A[, , 1] * m_R_effects 
  
  #### run the model
  for(t in 1:n_t){
    m_M[t + 1, ] <- m_M[t, ] %*% m_P             # estimate the state vector for the next cycle (t + 1)
    a_A[, , t + 1] <- diag(m_M[t, ]) %*% m_P     # estimate the transition dynamics at t + 1
    
    # element-wise - multiplication of array A with the rewards matrices to apply both state and transition rewards
    a_Y_costs[, , t + 1]   <- a_A[, , t + 1] * m_R_costs   
    a_Y_effects[, , t + 1] <- a_A[, , t + 1] * m_R_effects 
  }
  
  #### Aggregate outcomes ####
  v_costs   <- rowSums(t(colSums(a_Y_costs)))    # calculate the expected costs per cycle
  v_QALYs   <- rowSums(t(colSums(a_Y_effects)))  # calculate the expected QALYs per cycle
  TC        <- sum(v_costs)                      # calculate the total expected costs
  TE        <- sum(v_QALYs)                      # calculate the total expected QALYS
  v_results <- c(TC, TE)                         # combine the total expected costs and QALYs
  names(v_results) <- c("Costs", "Effect")       # name the vector

  #### memory
  # store all relevant objects and compute the memory it takes
  stored <- list(m_P, m_M, m_R_costs, m_R_effects, a_A, a_Y_costs, a_Y_effects)
  mem    <- object.size(stored)
  
  output <- list(results = v_results, memory = mem)
  return(output)
}

# test function
run_array_approach(n_states = 3, n_t = 60)

################ run simulations ####################

v_n_states <- seq(2, 62, 5)
# v_n_t <- seq(20, 1320, 100)  # (3 months to 1320 months (110 years)) 
v_n_t <- seq(12, 1320, 12)


# traditional
time_trace <- memory_trace <- c() # store time and memory

for (s in v_n_states) {
  for (t in v_n_t) {
    time_trace   <- c(time_trace, as.numeric(system.time(run_trace_approach(s, t))['elapsed']))
    memory_trace <- c(memory_trace, run_trace_approach(s, t)$memory)  
  }
}  

# array
time_array <- memory_array <- c() # store time and memory

for (s in v_n_states) {
  for (t in v_n_t) {
    time_array   <- c(time_array, as.numeric(system.time(run_array_approach(s, t))['elapsed']))
    memory_array <- c(memory_array, run_array_approach(s, t)$memory)  
  }
}  

##### simulation results #####

# cohort trace
time_mem_trace <- as.data.frame(cbind(rep("cohort trace", times = length(v_n_states) * length(v_n_t)), 
                                      rep(v_n_states, each = length(v_n_t)), 
                                      rep(v_n_t, times = length(v_n_states)),
                                      time_trace,
                                      memory_trace)) 
colnames(time_mem_trace) <- c('method', 'n_states', 'n_t', 'time', 'memory')

# array
time_mem_array <- as.data.frame(cbind(rep("array", times = length(v_n_states) * length(v_n_t)), 
                                      rep(v_n_states, each = length(v_n_t)), 
                                      rep(v_n_t, times = length(v_n_states)),
                                      time_array,
                                      memory_array)) 
colnames(time_mem_array) <- c('method', 'n_states', 'n_t', 'time', 'memory')

time_mem_both <- rbind(time_mem_trace, time_mem_array)

write.csv(time_mem_both, 'time_memory_comparisons_more_cycles.csv', row.names = FALSE)



