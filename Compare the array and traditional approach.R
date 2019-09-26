# Test if the array approach and traditional appraoch give similar results. 

# Run the array approach
source("Supplementary file_R code of the stylistic 3-state model.R")

# Store the array approach values in new variables
v_results_array_approach <- v_results
m_M_array_appraoch <- m_M

# Run the traditional approach 
source("Supplementary file_R code of the stylistic 3-state model_traditional approach.R")

# Store the results in a new variable
v_results_traditional_approach <- v_results

# Create a new trace matrix to sum the markov trace of the traditional approach. 
m_M_traditional_approach <- matrix(NA, 
                                    nrow = n_t + 1 , 
                                    ncol  = 3)
colnames(m_M_traditional_approach) <- colnames(m_M_array_approach) # name the columns

# Store the data in the right column
m_M_traditional_approach[, "H"] <- m_M[, "H"]
m_M_traditional_approach[, "S"] <- m_M[,"Stemp"] + m_M[, "S"] # sum the proportions of sick (Stemp + S)
m_M_traditional_appraoch[, "D"] <- m_M[,"Dtemp"] + m_M[, "D"] # sum the proportions of dead (Dtemp + D)

# Compare the two Markov traces 
round(m_M_traditional_appraoch, 10) == round(m_M_array_appraoch, 10)

# Compare the costs and effects
round(v_results_array_approach["Costs"], 5) == round(v_results_traditional_approach["Costs"], 5)
round(v_results_array_approach["Effect"], 5) == round(v_results_traditional_approach["Effect"], 5)


