#------------------------------------------------------#
#### Generate base-case parameters set              ####
#------------------------------------------------------#
f.generate_init_params <- function(){ # User defined
  v.params.init <- data.frame(
    # Costs
    c.H   = 2000,   # cost of remaining one cycle healthy 
    c.S1  = 4000,   # cost of remaining one cycle sick 
    c.S2  = 15000,  # cost of remaining one cycle sicker 
    c.D   = 0,      # cost of being dead (per cycle)
    #c.Trt = 12000,  # cost of treatment (per cycle 
    # State utilities
    u.H   = 1,      # utility when healthy 
    u.S1  = 0.75,   # utility when sick 
    u.S2  = 0.5,    # utility when sicker
    u.D   = 0,      # utility when healthy 
    #u.Trt = 0.95,  # utility when being treated
    # Transition probabilities (per cycle)
    p.HS1 = 0.15,   # probability to become sick when healthy
    p.S1H = 0.5,    # probability to become healthy when sick
    p.S1S2 = 0.105, # probability to become sicker when sick
    # Hazard rate ratios
    hr.S1  = 3,     # hazard rate ratio of death in S1 vs healthy
    hr.S2  = 10,    # hazard rateratio of death in S2 vs healthy 
    # Transition rewards
    du.HS1 = 0.01,  # one-time utility decrement when becoming sick
    ic.HS1 = 1000,  # one-time cost of becoming sick
    ic.D   = 2000   # one-time cost of dying
  )
  return(v.params.init)
}
