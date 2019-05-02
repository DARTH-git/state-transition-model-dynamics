# Check that the cumulative incidence is correct by checking the values from 
# the same model using 4 health states

##########

# Proof code
# Load the data from the 3-state model with the 4th temporary state 
load("Cohort_trace_4states.RData") # load the cohort trace 
load("Array_4states.RData") # load the array 


# Check if the 3-state model with 4 healht states gives the same results
round(m.M[, "H"], 12) == round(m.M_4states[, "H"] + m.M_4states[, "H2"], 12)
round(m.M[, "S"], 12) == round(m.M_4states[, "S"], 12)
round(m.M[, "D"], 12) == round(m.M_4states[, "D"], 12)
# Conclusion: YES, up to 12 decimals they are correct









n.atRisk <- a.A["H", "H", 1] * cumprod(p.HH)

n.atRisk <- m.M[1, "H"] * m.P["H", "H"]^(c(0, v.t[-length(v.t)]))  # calculate the total number at risk of getting sick for the first time (e.g. those that stayed out healthy since the start of the model). The first time point, cycle 0, everyone is at risk



# p.HS is equal for those that always stayed healthy and those that recovered from sick
n.HS <- a.A["H", "S", ]
n.H  <- colSums(a.A["H", , ])
p.HS <- n.HS / n.H  # get the transition probability at each time point from array A. Those that getting sick / total started healthy
n.newCases <- n.atRisk * p.HS  # new cases using the transition probability 
v.CI <- cumsum(n.newCases) / m.M[1, "H"] # calculate the cumulative incidence 



# Method 1 
for (k in 1:(n.t + 1)){ # for cycle 0 to the last cycle
  n.atRisk[k+1] <- n.atRisk[k] * m.P["H", "H"]
  n.atRisk[k] <- m.M[1, "H"] * (m.P["H", "H"]^(k - 1)) # proportion of the cohort at risk of getting sick for the first time # k - 1 : the first cycle they are all at risk
  
  p.HS[k - 1] <- a.A["H", "S", k] / sum(a.A["H", , k]) # get the transition probability at each time point from array A. Those that getting sick / total started healthy
  
  n.newCases[k + 1] <- n.atRisk[k] * p.HS[k]  # new cases using the transition probability 
  v.CI[k] <- sum(n.newCases[1:k])  / m.M[1, "H"] # calculate the cumulative incidence 
}

# Method 2 : using proportion 
for(k in 1:(n.t + 1)){
  n.atRisk[k] <- m.M[1, "H"] * (m.P["H", "H"]^(k - 1)) 
  p.atRisk[k + 1] <- n.atRisk[k] / m.M[k, "H"] # proportion of the cohort that is at risk of becomming sick for the first time from those that are healthy
  n.newCases2[k] <- a.A["H", "S", k] * p.atRisk[k]  # p.atRisk[k] # new of each time point
  v.CI_2[k] <- sum(n.newCases2[1:k])  / m.M[1, "H"] # calculate the cumulative incidence
}



# step 1: check if the calcualtion of the number at risk is correct 
round(n.atRisk, 12) == round(m.M_4states[, "H"], 12) # they are equal up to 15 decimals
# Conclusion: calculation of the number at risk is correct.

# step 2: check if the number of healthy individuals is correct 
round(m.M[, "H"], 12) == round(m.M_4states[, "H"] +  m.M_4states[, "H2"], 12) 
# Conclusion: calculation of those that are healthy is correct up to 14 decimals

# step 3: check if the ratio of at risk/healthy is correct 
p.atRisk_4states <- m.M_4states[, "H"] / (m.M_4states[, "H"] +  m.M_4states[, "H2"])
round(p.atRisk, 12) == round(p.atRisk_4states, 12) 


# Conclusion : they are correct, but not with 15 decimals anymore, but with 14
# This is still good enough.

# step 4: check if the total number of individuals that transition to sick is correct
round(a.A["H", "S", ], 12)  == round(a.A_4states["H", "S", ] + a.A_4states["H2", "S", ], 12)
# total number of individuals that move from health to sick 
# conclussion = correct 

# step 5: check if the vector with the cumulative incidence is equal to the
round(v.CI, 12) == round(cumsum(a.A_4states["H", "S", ]), 12) # 


