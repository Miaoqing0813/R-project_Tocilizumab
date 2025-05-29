#######################
# Semi-Markov model #
######################

# Markov States: 
n_states <- 4  
state_names <- c("Initial remission", "Relapse", "Subsequent remission", "Dead")

# Time-dependent transition probability matrices (TPM) creation function
get_tpm <- function(arm, cycle) {
  tpm <- matrix(0, nrow = n_states, ncol = n_states,
                dimnames = list(from = state_names, to = state_names))
  
  p_allcause_mort <- v.allcause.mort.monthly[cycle]
  
  if (arm == "TCZ_plus_GC") {
    p_first_flare <- if (cycle <= 26) p_first_flare_TCZ_0_26 else p_first_flare_TCZ_27_end
    p_flare <- exp(log_flare_TCZ)
  } else {
    p_first_flare <- if (cycle <= 26) p_first_flare_GC_0_26 else p_first_flare_GC_27_end
    p_flare <- exp(log_flare_GC)
  }
  
  # From Initial remission
  tpm["Initial remission", "Relapse"] <- p_first_flare
  tpm["Initial remission", "Dead"] <- p_allcause_mort
  tpm["Initial remission", "Initial remission"] <- 1 - p_first_flare - p_allcause_mort
  
  # From Relapse: move directly to remission if alive
  tpm["Relapse", "Subsequent remission"] <- 1 - p_allcause_mort
  tpm["Relapse", "Dead"] <- p_allcause_mort
  tpm["Relapse", "Relapse"] <- 0
  tpm["Relapse", "Initial remission"] <- 0
  
  # From Subsequent remission
  tpm["Subsequent remission", "Relapse"] <- p_flare
  tpm["Subsequent remission", "Dead"] <- p_allcause_mort
  tpm["Subsequent remission", "Subsequent remission"] <- 1 - p_flare - p_allcause_mort
  
  # Absorbing state
  tpm["Dead", "Dead"] <- 1
  
  return(tpm)
}

# state populations 
pop <- array(data = NA,
             dim = c(n_states, n_cycles, n_arms),
             dimnames = list(state = state_names,
                             cycle = NULL,
                             arms = arms))

# Assume all patients starts with initial remission: newly diagnosed patients
pop["Initial remission", cycle = 1, ] <- 1
pop[-1, cycle = 1, ] <- 0

# New patients: useful if Y1 cost is different from Y2 onwards
trans <- array(data = NA,
               dim = c(n_states, n_cycles, n_arms),
               dimnames = list(state = state_names,
                               cycle = NULL,
                               arms = arms))

# Sum costs and QALYs for each cycle at a time for each diagnostic outcome 
cycle_empty_array <- array(NA,
                           dim = c(n_arms, n_cycles),
                           dimnames = list(arms = arms,
                           cycle = NULL))

cycle_costs <- cycle_qalys <- cycle_empty_array

