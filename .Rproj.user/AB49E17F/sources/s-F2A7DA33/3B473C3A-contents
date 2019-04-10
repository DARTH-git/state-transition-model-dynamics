#------------------------------------------------------#
#### Generate the transition probability matrix     ####
#------------------------------------------------------#
f.create_transition_prob_matrix <- function(df.params){# User defined
  ### Description:
  ##    This function constructs the transition probability matrix m.P
  ### Arguments:  
  ##    df.params: dateframe with the model parameters 
  ### Returns
  ##     m.P: the transition probability matrix 
  with(as.list(df.params), {
    # Initialize matrix
    m.P <- matrix(NA, 
                  nrow = n.states, ncol = n.states, 
                  dimnames = list(v.n, v.n))
    
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
    
    return(m.P) # return the transition probability matrix 
  }
  )
}

#------------------------------------------------------#
#### Generate the reward matrix for costs          ####
#------------------------------------------------------#
f.create_transition_reward_matrix_costs <- function(df.params){ # User defined
  ### Description:
  ##    This function generates the reward matrix for costs
  ### Arguments:  
  #     df.params: dataframe with the model parameters 
  ### Returns
  #     m.R: the reward matrix 
  with(as.list(df.params), {
    
    # initialize matrix
    m.R <- matrix(NA, 
                  nrow = n.states, ncol = n.states, 
                  dimnames = list(v.n, v.n))
    
    # Fill in matrix
    # From Healthy
    m.R["H", "H"]   <- c.H
    m.R["H", "S"]   <- c.H 
    m.R["H", "D"]   <- c.H + ic.D
    # From Sick
    m.R["S", "H"]  <- c.S
    m.R["S", "S"]  <- c.S 
    m.R["S", "D"]  <- c.S + ic.D
    # From Death
    m.R["D", "H"]   <- c.D
    m.R["D", "S"]   <- c.D
    m.R["D", "D"]   <- c.D 
    
    return(m.R) # return the reward matrix 
  }
  )
}

#------------------------------------------------------#
#### Generate the reward matrix for effects         ####
#------------------------------------------------------#
f.create_transition_reward_matrix_effects <- function(df.params){ # User defined
  ### Description:
  ##    This function generates the reward matrix for effects
  ### Arguments:   
  #     v.params: vector of model parameters 
  ### Returns
  #     m.R: the reward matrix 
  with(as.list(df.params), {
    
    # initialize matrix
    m.R <- matrix(NA, 
                  nrow = n.states, ncol = n.states, 
                  dimnames = list(v.n, v.n))
    
    # Fill in matrix
    # From Healthy
    m.R["H", "H"]   <- u.H
    m.R["H", "S"]   <- u.H - du.HS
    m.R["H", "D"]   <- u.H 
    # From Sick
    m.R["S", "H"]   <- u.S
    m.R["S", "S"]   <- u.S
    m.R["S", "D"]   <- u.S 
    # From Death
    m.R["D", "H"]   <- u.D
    m.R["D", "S"]   <- u.D
    m.R["D", "D"]   <- u.D 
    
    return(m.R) # return the reward matrix 
  }
  )
}