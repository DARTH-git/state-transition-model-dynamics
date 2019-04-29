###############################################################################
### A brief report -  ## 2019 ##
################################################################################
# This code forms the basis for the brief report: 
# An alternative representation of state transition model dynamics.
# Authors: Eline Krijkamp^{co}, Fernando Alarid-Escudero^{co}, Eva A. Enns, 
# Myriam G.M. Hunink, Petros Pechlivanoglou, Hawre Jalal.

# Please cite the article when using this code
################################################################################
# Demonstrate the array appraoch using the Sick-Sicker model with age dependent
# transition probabilities but without performing a cost-effectiveness analysis.
################################################################################
# To program this tutorial we made use of 
# R version 3.5.1 (2018-7-02)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14
# RStudio: Version 1.1.456 2009-2018 RStudio, Inc
################################################################################
rm(list = ls())  # remove any variables in R's memory 

###01 Initial setup 
#### 01.1 Load packages and functions ####
library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(scales)   # for dollar signs and commas
library(tensorA)  # for tensor calculations 

#### 01.1.2 Load functions ####
source("Sick-Sicker case example/functions/01_model-inputs_functions.R")
source("Sick-Sicker case example/functions/02_simulation-model_functions.R")

#### 01.2 External parameters ####
#### 01.2.1 General setup ####
n.age.init  <- 25  # age of starting cohort
n.t         <- 75  # time horizon, number of cycles
v.age.names <- n.age.init:(n.age.init + n.t - 1) # vector with age names
v.n <- c("H", "S1", "S2", "D") # vector with the 4 health states of the model:
# Healthy (H), Sick (S1), Sicker (S2), Dead (D)
n.states <- length(v.n) # number of health states 
d.c <- 0.03 # discount rate for costs 
d.e <- 0.03 # discount rate for QALYs
v.dwc <- 1 / ((1 + d.e) ^ (0:(n.t))) # vector with discount weights for costs
v.dwe <- 1 / ((1 + d.c) ^ (0:(n.t))) # vector with discount weights for QALYs
v.s.init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector

#### 01.2.2 All-cause age-, sex- and race- (ASR) specific mortality ####
df.r.asr <- read.csv("Sick-Sicker case example/data/01_all-cause-mortality-USA-2015.csv")
v.r.asr  <- df.r.asr %>%
  dplyr::select(Total) %>%
  as.matrix()              # vector with mortality rates

#### 01.2.3 Generate initial set of base-case external parameters ####
v.params.init <- f.generate_init_params()
## Create name of parameters
v.names.params <- names(v.params.init)


### 02 Define and initialize matrices and vectors ####
#### 02.1 Transition probability matrix ####
#### Equation 1 #### 
f.create_transition_prob_matrix(v.params = v.params.init, t = 1) # the m.P transition probability matrix at the first cycle

#### 02.2 Initial state vector ####
# The cohort start in the Healthy health state
s0 <- c(H = 1, S1 = 0, S2 = 0, D = 0)
s0

#### 02.3 Cohort trace  
## Create the Markov trace matrix M capturing the proportion of the cohort in each state at each cycle
# Initialize cohort trace
m.M <- matrix(0, 
              nrow = (n.t + 1), ncol = n.states, 
              dimnames = list(0:n.t, v.n))

m.M[1, ] <- s0 # store the initial state vector


#### 03 Matrix Approach    ####
for(t in 1:n.t){  # loop through the number of cycles
  m.P <- f.create_transition_prob_matrix(v.params = v.params.init, t = t) # create the transition probability matrix for the current cycle
#### Equation 2   #### 
  m.M[t + 1, ] <- m.M[t, ] %*% m.P  # estimate the state vector for the next cycle (t + 1)
}

head(round(m.M, 3)) # show the first six lines of the Markov trace

#### 04 Array Approach     ####
a.A <- array(0, dim = c(n.states, n.states, n.t + 1),
             dimnames = list(v.n, v.n, 0:n.t)) # Initialize array
diag(a.A[, , 1]) <- s0 # store the initial state vector in the diagnal of A
#### Equation 3  #### 
a.A[, , 1]

# run the model 
for(t in 1:n.t){                     # loop through the number of cycles
  m.P <- f.create_transition_prob_matrix(v.params = v.params.init, t = t) # create the transition probability matrix for the current cycle
#### Equation 4    #### 
  a.A[, , t + 1] <- colSums(a.A[, , t]) * m.P  # fill array A for t + 1 
}

#### Equation 5    #### 
a.A[, , 2:3] # shown for two cycles

#### Equation 7    #### 
# calculating M from A 
m.M_A <- t(colSums(a.A))   # sum over the colums of A and transpose 
m.M == m.M_A # check if they are exactly the same


#### 05 Apply state and transtion rewards #### 
#### 05.1 Create reward matrices for both costs and effects #### 
m.R_costs  <- f.create_transition_reward_matrix_costs(v.params = v.params.init)
m.R_effects <- f.create_transition_reward_matrix_effects(v.params = v.params.init)
#### Equation 8   #### 
m.R_costs 
m.R_effects 


#### 05.2 Expected QALYs and Costs per cycle for each strategy ####
#### Equation 9 ####
a.Y_costs <- a.Y_effects <- array(0, dim = c(n.states, n.states, n.t + 1),
             dimnames = list(v.n, v.n, 0:n.t))

for(t in 1:n.t){ 
a.Y_costs[, , t]   <- a.A[, , t] * m.R_costs
a.Y_effects[, , t] <- a.A[, , t] * m.R_effects
}

## Vector of expected costs per cycle
#v.cost_UC  <- rowSums(t(colSums(to.tensor(a.A) * to.tensor(m.R_costs))))
## Vector of expected QALYs per cycle
#v.qaly_UC  <- rowSums(t(colSums(to.tensor(a.A) * to.tensor(m.R_effect))))

#### Equation 10 ####
v.Costs <- rowSums(t(colSums(a.Y_costs)))
v.QALYs <- rowSums(t(colSums(a.Y_effects)))

TC <- t(v.Costs) %*% v.dwc
TE <- t(v.QALYs) %*% v.dwe

v.Results <- c(TC, TE)
names(v.Results) <- c("Costs", "Effect")
v.Results # print the results


#### 06 Plot cohort trace  ####
ggplot(melt(m.M), aes(x = Var1, y = value, color = Var2)) +
  geom_line(size = 1.3) +
  scale_color_discrete(l = 50, name = "Health state", h = c(45, 365)) +
  xlab("Cycle") +
  ylab("Proportion of the cohort") +
  theme_bw(base_size = 16) +
  theme()

