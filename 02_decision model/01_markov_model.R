#######################
# Semi-Markov model #
######################

# Markov States: 
n_states <- 5  
state_names <- c("Active GCA", "Initial remission", "Relapse", "Subsequent remission", "Dead")

# Time-dependent transition probability matrices (TPM) creation function
get_tpm <- function(arm, cycle) {
  
  # Ensure cycle does not exceed the length of input data
  cycle <- min(cycle, 360) # monthly cycle for 30 years
  
  # Get monthly all-cause mortality rate with steroid hazard ratio adjustment
  p_allcause_mort <- get_steroid_ae_probs(cycle = cycle)$mortality
  
  # Extract cycle-specific flare probabilities
  if (arm == 1) {
    p_first_flare <- trans_probs_first_flare_treat_1yr_weibull$trans_prob[cycle]
    p_flare <- probs_subseq_flare_treat
  } else {
    p_first_flare <- trans_probs_first_flare_control_1yr_weibull$trans_prob[cycle]
    p_flare <- probs_subseq_flare_control
  }
  
  # Initialize TPM
  tpm <- matrix(0, nrow = n_states, ncol = n_states,
                dimnames = list(from = state_names, to = state_names))
  
  # From Active GCA
  if (cycle < 13) {
    # Within first year
    tpm["Active GCA", "Relapse"] <- p_first_flare
    tpm["Active GCA", "Dead"] <- p_allcause_mort
    tpm["Active GCA", "Active GCA"] <- 1 - p_first_flare - p_allcause_mort
  } else {
    # After first year, move to remission if still active
    tpm["Active GCA", "Initial remission"] <- 1 - p_first_flare - p_allcause_mort
    tpm["Active GCA", "Relapse"] <- p_first_flare
    tpm["Active GCA", "Dead"] <- p_allcause_mort
  }
  
  # From Initial remission
  if (cycle >= 13) {
    tpm["Initial remission", "Relapse"] <- p_first_flare
    tpm["Initial remission", "Dead"] <- p_allcause_mort
    tpm["Initial remission", "Initial remission"] <- 1 - p_first_flare - p_allcause_mort
  } else {
    # Before cycle 53 nobody should be in remission
    tpm["Initial remission", "Initial remission"] <- 1
  }
  
  # From Relapse
  tpm["Relapse", "Subsequent remission"] <- 1 - p_allcause_mort
  tpm["Relapse", "Dead"] <- p_allcause_mort
  
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

# Assume all patients starts with Active GCA: newly diagnosed patients
pop["Active GCA", cycle = 1, ] <- 1
pop[-1, cycle = 1, ] <- 0

# New patients: useful if Y1 cost is different from Y2 onwards
trans <- array(data = NA,
               dim = c(n_states, n_cycles, n_arms),
               dimnames = list(state = state_names,
                               cycle = NULL,
                               arms = arms))

# Model GCA-related complications and steroid-related adverse events as events rather than distinct health states
event_log <- list(
  vl = array(0, dim = c(n_arms, n_cycles)),
  stroke = array(0, dim = c(n_arms, n_cycles)),
  mi = array(0, dim = c(n_arms, n_cycles)),
  hf = array(0, dim = c(n_arms, n_cycles)),
  fractures = array(0, dim = c(n_arms, n_cycles)),
  diabetes = array(0, dim = c(n_arms, n_cycles)), 
  infection = array(0, dim = c(n_arms, n_cycles))
)

# Track cumulative populations for all conditions
cumulative_populations <- list(
  vl = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  stroke = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  mi = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  hf = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  fractures = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  diabetes = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  infection = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL))
)

# Sum costs and QALYs for each cycle at a time for each diagnostic outcome 
cycle_empty_array <- array(0,
                           dim = c(n_arms, n_cycles),
                           dimnames = list(arms = arms,
                           cycle = NULL))

cycle_costs <- cycle_qalys <- cycle_empty_array
event_costs <- event_utility <- cycle_empty_array

# Running the model 

## Update pop -- number of people at each state per cycle
for (i in 1:n_arms) {
  for (j in 2:n_cycles) {
    tpm <- get_tpm(arm = i, cycle = j - 1)
    pop[, cycle = j, arm = i] <- pop[, cycle = j - 1, arm = i] %*% tpm
    
    # Loop over states and apply event probabilities
    for (state in state_names[state_names != "Dead"]) {
      
      state_pop <- pop[state, j, i]
      
      # Define cycle time in years
      t_years <- j * cycle_length / 12
      current_year <- floor(t_years) + 1  # Start from year 1, not 0
      
      # GCA-related complications: only in Relapse
      if (state == "Relapse") {
        # Vision loss - only new incidence within year 1 (cycle 13)
        if (j <= 13) {
          n_vl_incident <- state_pop * p_vl
          event_log$vl[i, j] <- n_vl_incident
          
          # Add new vision loss cases to cumulative population
          if (j == 1) {
            cumulative_populations$vl[i, j] <- n_vl_incident
          } else {
            cumulative_populations$vl[i, j] <- cumulative_populations$vl[i, j-1] + n_vl_incident
          }
        } else {
          # After year 1, carry forward the vision loss population (no new cases)
          cumulative_populations$vl[i, j] <- cumulative_populations$vl[i, j-1]
        }
      } else {
        # For non-Relapse states, carry forward the vision loss population
        if (j > 1) {
          cumulative_populations$vl[i, j] <- cumulative_populations$vl[i, j-1]
        }
      }
      
      # Steroid-related events (all non-death states) 
      steroid_ae_probs <- get_steroid_ae_probs(cycle = j)
      
      for (event in names(steroid_ae_probs)) {
        if (event != "mortality") {  # Skip mortality as it's handled separately
          prob_event <- steroid_ae_probs[[event]]
          n_event <- state_pop * prob_event
          event_log[[event]][i, j] <- event_log[[event]][i, j] + n_event
          
          # Update cumulative population for this event
          if (j == 1) {
            cumulative_populations[[event]][i, j] <- n_event
          } else {
            cumulative_populations[[event]][i, j] <- cumulative_populations[[event]][i, j-1] + n_event
          }
        }
      }
    }
    
    # After processing all states, apply costs and utility impacts for all permanent conditions
    # Vision loss costs and utility impacts
    if (cumulative_populations$vl[i, j] > 0) {
      vl_cost <- ifelse(current_year == 1, c_vl_y1, c_vl_y2on)
      event_costs[i, j] <- event_costs[i, j] + apply_discount(cumulative_populations$vl[i, j] * vl_cost, disc_rate, t_years)
      event_utility[i, j] <- event_utility[i, j] + apply_discount(cumulative_populations$vl[i, j] * q_vl, disc_rate, t_years)
    }
    
    # Steroid-related adverse events costs and utility impacts
    steroid_events <- c("stroke", "mi", "hf", "fractures", "diabetes", "infection")
    for (event in steroid_events) {
      if (cumulative_populations[[event]][i, j] > 0) {
        # Event-specific cost and utility loss
        if (event == "stroke") {
          cost <- ifelse(current_year == 1, c_stroke_y1, c_stroke_y2on)
          qloss <- ifelse(current_year == 1, q_stroke_y1, q_stroke_y2on)
        } else if (event == "mi") {
          cost <- ifelse(current_year == 1, c_mi_y1, c_mi_y2on)
          qloss <- ifelse(current_year == 1, q_mi_y1, q_mi_y2on)
        } else if (event == "hf") {
          cost <- ifelse(current_year == 1, c_hf_y1, c_hf_y2on)
          qloss <- ifelse(current_year == 1, q_hf_y1, q_hf_y2on)
        } else if (event == "fractures") {
          cost <- ifelse(current_year == 1, c_fracture_y1, c_fracture_y2on)
          qloss <- ifelse(current_year == 1, q_fracture_y1, q_fracture_y2on)
        } else if (event == "diabetes") {
          cost <- c_dm
          qloss <- q_dm
        } else if (event == "infection") {
          cost <- c_infection
          qloss <- q_infection
        } else {
          cost <- 0
          qloss <- 0
        }
        
        event_costs[i, j] <- event_costs[i, j] + apply_discount(cumulative_populations[[event]][i, j] * cost, disc_rate, t_years)
        event_utility[i, j] <- event_utility[i, j] + apply_discount(cumulative_populations[[event]][i, j] * qloss, disc_rate, t_years)
      }
    }
  }
}
