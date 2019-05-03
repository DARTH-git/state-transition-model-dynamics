# Check that the cumulative incidence is correct by checking the values from 
# the same model using 4 health states

##########
rm(list = ls())  # remove any variables in R's memory 

# Load the data
## 3-state model
load("output/Cohort_trace.RData") # load the cohort trace 
load("output/Array.RData") # load the array
load("output/CumulativeIncidence.RData") # load the cumulative incidence 
load("output/n.atRisk.RData") # load the number at risk
load("output/n.newCases.RData") # load the number of new cases 

## 3-state model with the 4th temporary state 
load("output/Cohort_trace_4states.RData") # load the cohort trace 
load("output/Array_4states.RData") # load the array 

# Step 1: check if the 3-state model with 4 health states gives the same results
round(m.M[, "H"], 12) == round(m.M_4states[, "H"] + m.M_4states[, "H2"], 12)
round(m.M[, "S"], 12) == round(m.M_4states[, "S"], 12)
round(m.M[, "D"], 12) == round(m.M_4states[, "D"], 12)
# Conclusion: Yes, identical up to 12 decimals 

# Step 2: check if the calculation of the number at risk is correct 
round(n.atRisk, 12) == round(m.M_4states[, "H"], 12) 
# Conclusion: Yes, identical up to 12 decimals 

# Step 3: check if the total number of individuals that transition from healthy to sick is correct
round(a.A["H", "S", ], 12)  == round(a.A_4states["H", "S", ] + a.A_4states["H2", "S", ], 12)
# total number of individuals that move from health to sick 
# Conclusion: Yes, identical up to 12 decimals 

# Step 4: check if the number of new cases is correct
round(n.newCases, 12) == round(a.A_4states["H", "S", ], 12) # 
# Conclusion: Yes, identical up to 12 decimals 

# Step 5: check if the vector with the cumulative incidence is equal to the the cumulative incidence based on the 4 state model
round(v.CI, 12) == round(cumsum(a.A_4states["H", "S", ]), 12) # 
# Conclusion: Yes, identical up to 12 decimals 

################################################################################



