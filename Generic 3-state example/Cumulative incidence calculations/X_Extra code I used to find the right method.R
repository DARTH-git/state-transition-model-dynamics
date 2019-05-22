







n.atRisk <- a.A["H", "H", 1] * cumprod(p.HH)

n.atRisk <- m.M[1, "H"] * m.P["H", "H"]^(c(0, v.t[-length(v.t)]))  # calculate the total number at risk of getting sick for the first time (e.g. those that stayed out healthy since the start of the model). The first time point, cycle 0, everyone is at risk



# p.HS is equal for those that always stayed healthy and those that recovered from sick
n.HS <- a.A["H", "S", ]
n.H  <- colSums(a.A["H", , ])
p.HS <- n.HS / n.H  # get the transition probability at each time point from array A. Those that getting sick / total started healthy
n.newCases <- n.atRisk * p.HS  # new cases using the transition probability 
v.CI <- cumsum(n.newCases) / m.M[1, "H"] # calculate the cumulative incidence 


# step 3: check if the ratio of at risk/healthy is correct 
p.atRisk_4states <- m.M_4states[, "H"] / (m.M_4states[, "H"] +  m.M_4states[, "H2"])
round(p.atRisk, 12) == round(p.atRisk_4states, 12) 

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

