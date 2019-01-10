###############################################################################
### A brief report -  ## 2019 ##
################################################################################
# This code forms the basis for the brief report: 
# 'State transition model dynamics TITLE' 
# Authors: 
# Please cite the article when using this code
#
# To program this tutorial we made use of 
# R version 3.5.1 (2018-7-02)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14
# RStudio: Version 1.1.456 2009-2018 RStudio, Inc

################################################################################
################# Code of Appendix A ###########################################
################################################################################
# Constructs a time-homogeneous implementation of the Sick-Sicker model #

################################# Initial setup ################################
#rm(list = ls())  # remove any variables in R's memory 
library(dplyr)    # to manipulate data
library(reshape2) # to transform data
library(ggplot2)  # for nice looking plots
library(scales)   # for dollar signs and commas

#### 02 Load Functions ####
# no functions required

#### 03 Input Model Parameters ####
## Strategy names
v.names.str <- c("No Treatment", "Treatment")  
## Number of strategies
n.str <- length(v.names.str)
## Markov model parameters
age     <- 25                  # age at baseline
max.age <- 55                  # maximum age of follow up
n.t <- max.age - age           # time horizon, number of cycles
v.n <- c("H", "S1", "S2", "D") # the 4 states of the model: 
# Healthy (H), Sick (S1), Sicker (S2), Dead (D)
n.s <- length(v.n)             # number of health states 

## Transition probabilities (per cycle) and hazard ratios
p.HD    <- 0.002 # constant probability of dying when healthy (all-cause mortality)
p.HS1   <- 0.15  # probability to become sick when healthy
p.S1H   <- 0.5   # probability to become healthy when sick
p.S1S2  <- 0.105 # probability to become sicker when sick
hr.S1   <- 3     # hazard ratio of death in S1 vs healthy
hr.S2   <- 10    # hazard ratio of death in S2 vs healthy 

## Cost and utility inputs 
# State rewards
c.H   <- 2000  # cost of remaining one cycle healthy 
c.S1  <- 4000  # cost of remaining one cycle sick 
c.S2  <- 15000 # cost of remaining one cycle sicker 
c.D   <- 0     # cost of being dead (per cycle)
c.Trt <- 12000 # cost of treatment (per cycle 

u.H   <- 1     # utility when healthy 
u.S1  <- 0.75  # utility when sick 
u.S2  <- 0.5   # utility when sicker
u.D   <- 0     # utility when healthy 
u.Trt <- 0.95  # utility when being treated

# Transition rewards
du.HS1 <- 0.01  # disutility when transitioning from H to S1
ic.HS1 <- 1000  # increase in cost when transitioning from H to S1
ic.D   <- 1000  # increase in cost when dying

## Discounting factor
d.c <- d.e <- 0.03 # equal discount of costs and QALYs by 3%
# Discount weight vector (equal discounting is assumed for costs and effects)
v.dwc <- 1 / (1 + d.c) ^ (0:n.t) # calculate discount weights for costs for each cycle based on discount rate d.r
v.dwe <- 1 / (1 + d.e) ^ (0:n.t) # calculate discount weights for effectiveness for each cycle based on discount rate d.r

# Transition probability to dead from S1 and S2
p.S1D  <- 1 - exp(log(1 - p.HD) * hr.S1) # probability to die in S1
p.S2D  <- 1 - exp(log(1 - p.HD) * hr.S2) # probability to die in S2

#### 04 Define and initialize matrices and vectors ####
#### 04.1 Transition probability MATRIX ####
# Initialize matrix
m.P <- matrix(0, 
              nrow = n.s, ncol = n.s, 
              dimnames = list(v.n, v.n))
# Fill in matrix
# From H
m.P["H", "H"]  <- 1 - (p.HS1 + p.HD)
m.P["H", "S1"] <- p.HS1
m.P["H", "D"]  <- p.HD
# From S1
m.P["S1", "H"]  <- p.S1H
m.P["S1", "S1"] <- 1 - (p.S1H + p.S1S2 + p.S1D)
m.P["S1", "S2"] <- p.S1S2
m.P["S1", "D"]  <- p.S1D
# From S2
m.P["S2", "S2"] <- 1 - p.S2D
m.P["S2", "D"]  <- p.S2D
# From D
m.P["D", "D"] <- 1

#### 04.2 Initial state vector ####
# All starting healthy
s0 <- c(H = 1, S1 = 0, S2 = 0, D = 0)
s0

#### 04.3 Cohort trace ####
## Create the markov trace matrix M capturing the proportion of the cohort in each state at each cycle
# Initialize cohort trace
m.M <- matrix(0, 
              nrow = (n.t + 1), ncol = n.s, 
              dimnames = list(0:n.t, v.n))
m.M[1, ] <- s0 # store the initial state vector

################################################################################
############################## Matrix approach #################################
################################################################################
#### 05 Run Markov model ####
for(t in 1:n.t){                    # loop through the number of cycles
  m.M[t + 1, ] <- m.M[t, ] %*% m.P  # estimate the state vector for the next cycle (t + 1)
}
head(round(m.M, 3))

################################################################################
############################## Array approach ##################################
################################################################################
# Initialize array
a.A <- array(0, dim = c(n.s, n.s, n.t + 1),
             dimnames = list(v.n, v.n, 0:n.t))

#### 05.1 Run Markov model ####
# as described in the paper
t05.1 <- Sys.time()
diag(a.A[, , 1]) <- s0 # store the initial state vector in the diagnal of A

# run the model 
for(t in 1:n.t){                     # loop through the number of cycles
  a.A[, , t + 1] <- colSums(a.A[, ,t]) * m.P  # fill array A for t + 1 
}

# calculating M from A 
m.M_A <- t(colSums(a.A)) # sum over the colums of A and transpose 

t05.1 = Sys.time() - t05.1

#### 05.2 Run Markov model ####
# creating both the matrix and array at the same time 
t05.2 <- Sys.time()
for(t in 1:n.t){                     # loop through the number of cycles
  m.M[t + 1, ]   <- m.M[t, ] %*% m.P # estimate the state vector for cycle t + 1
  a.A[, , t + 1] <- m.M[t, ]  * m.P  # fill array A for t + 1 
}
t05.2 <- Sys.time() - t05.2

m.M_A == m.M # do the two approached give the same results?
t05.1 < t05.2 # is the first approach faster?


#### 06 Compute and Plot  ####
#### 06.1 Cohort trace #####
ggplot(melt(m.M), aes(x = Var1, y = value, color = Var2)) +
  geom_line(size = 1.3) +
  scale_color_discrete(l = 50, name = "Health state", h = c(45, 365)) +
  xlab("Cycle") +
  ylab("Proportion of the cohort") +
  theme_bw(base_size = 16) +
  theme()


#### 07 Compute Cost-Effectiveness Outcomes ####
#### 07.1 State rewards for each strategy ####
## Vector of state utilities under usual care
m.u.UC <- matrix(0, 
                 nrow = n.s, ncol = n.s, 
                 dimnames = list(v.n, v.n))

# From H
m.u.UC["H", "H"]    <- u.H
m.u.UC["H", "S1"]   <- u.S1 - du.HS1
m.u.UC["H", "D"]    <- u.D
# From S1
m.u.UC["S1", "H"]   <- u.H
m.u.UC["S1", "S1"]  <- u.S1 
m.u.UC["S1", "S2"]  <- u.S2
m.u.UC["S1", "D"]   <- u.D
# From S2
m.u.UC["S2", "S1"]  <- u.S1
m.u.UC["S2", "S2"]  <- u.S2
m.u.UC["S2", "D"]   <- u.D

## Matrix of state costs under usual care
m.c.UC <- matrix(0, 
                 nrow = n.s, ncol = n.s, 
                 dimnames = list(v.n, v.n))
# Fill in matrix
# From H
m.c.UC["H", "H"]  <- c.H 
m.c.UC["H", "S1"] <- c.S1 + ic.HS1
m.c.UC["H", "D"]  <- ic.D + c.D
# From S1
m.c.UC["S1", "H"]  <- c.H
m.c.UC["S1", "S1"] <- c.S1
m.c.UC["S1", "S2"] <- c.S2
m.c.UC["S1", "D"]  <- ic.D
# From S2
m.c.UC["S2", "S2"] <- c.S2
m.c.UC["S2", "D"]  <- ic.D
# From D
m.c.UC["D", "D"] <- c.D

## Vector of state costs under new treatment
m.c.Tr = m.c.UC # copy the results
# replace the cost values influences by treatment costs 
# From S1
m.c.Tr["S1", "S1"] <- c.S1 + c.Trt
m.c.Tr["S1", "S2"] <- c.S2 + c.Trt
# From S2
m.c.Tr["S2", "S2"] <- c.S2 + c.Trt


#### 07.2 Expected QALYs and Costs per cycle for each strategy ####
#### Expected QALYs and Costs per cycle ####
## Vector of qalys
v.qaly.UC <- rowSums(t(colSums(to.tensor(a.A) * to.tensor(m.u.UC))))
## Vector of costs
v.cost.UC <- rowSums(t(colSums(to.tensor(a.A) * to.tensor(m.c.UC))))


## Vector of qalys
v.qaly.Tr <- rowSums(t(colSums(to.tensor(a.A) * to.tensor(m.u.Tr))))
## Vector of costs
v.cost.Tr <- rowSums(t(colSums(to.tensor(a.A) * to.tensor(m.c.Tr))))

#### Discounted Total expected QALYs and Costs ####
te.UC <- v.qaly.UC %*% v.dwe  # total (discounted) QALY 
tc.UC <- v.cost.UC %*% v.dwc  # total (discounted) cost 

te.Tr <- v.qaly.Tr %*% v.dwe  # total (discounted) QALY 
tc.Tr <- v.cost.Tr %*% v.dwc  # total (discounted) cost 


#### Cost-effectiveness analysis ####
### Vector of costs
v.cost <- c(tc.UC, tc.Tr)
### Vector of effectiveness
v.qaly <- c(te.UC, te.Tr)

### Incremental outcomes
delta.C <- v.cost[2] - v.cost[1]             # calculate incremental costs
delta.E <- v.qaly[2] - v.qaly[1]             # calculate incremental QALYs
ICER    <- delta.C / delta.E                 # calculate the ICER
results <- c(delta.C, delta.E, ICER)         # store the values in a new variable

# Create full incremental cost-effectiveness analysis table
table_cstm <- data.frame(
  Costs = dollar(v.cost),
  QALYs = round(v.qaly, 3),
  `Incremental Costs` = c("", dollar(delta.C)),
  `Incremental QALYS` = c("", round(delta.E, 3)),
  ICER = c("", dollar(ICER))
)
rownames(table_cstm) <- v.names.str

table_cstm  # print the table 
