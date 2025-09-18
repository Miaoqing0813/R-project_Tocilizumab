#######################
# Semi-Markov model with Tunnel States #
######################

# For "Subsequent remission", we need tunnel states to track time-in-state
# Create tunnel states for subsequent remission (up to maximum duration needed)
max_tunnel_states <- steroid_effect_duration_months  # Track up to steroid effect duration
subseq_tunnel_names <- paste0("Subsequent remission_", 1:max_tunnel_states)

# Expanded state names including tunnel states
state_names <- c("Initial remission", "Relapse", subseq_tunnel_names, "Dead")
n_states <- length(state_names)

# Time-dependent transition probability matrices (TPM) creation function
get_tpm <- function(arm, cycle) {
  
  # Ensure cycle does not exceed the length of input data
  cycle <- min(cycle, n_cycles) 
  
  # Extract cycle-specific flare probabilities
  if (arm == 1) {
    p_first_flare <- trans_probs_first_flare_treat_weibull$trans_prob[cycle]
    p_flare <- probs_subseq_flare_treat
  } else {
    p_first_flare <- trans_probs_first_flare_control_weibull$trans_prob[cycle]
    p_flare <- probs_subseq_flare_control
  }
  
  # Initialize TPM
  tpm <- matrix(0, nrow = n_states, ncol = n_states,
                dimnames = list(from = state_names, to = state_names))
  
  # Get mortality for Initial remission (time_in_state = cycle number)
  mort_init <- get_steroid_ae_probs_by_state(
    cycle = cycle,
    state = "Initial remission",
    arm   = arm,
    duration = steroid_effect_duration_months,
    time_in_state = cycle
  )$mortality
  
  # From Initial remission
  tpm["Initial remission", "Relapse"]           <- p_first_flare
  tpm["Initial remission", "Dead"]              <- mort_init
  tpm["Initial remission", "Initial remission"] <- 1 - p_first_flare - mort_init
  
  # Get mortality for Relapse
  mort_relapse <- get_steroid_ae_probs_by_state(
    cycle = cycle,
    state = "Relapse",
    arm   = arm,
    duration = steroid_effect_duration_months,
    time_in_state = 1   # relapse always only one cycle
  )$mortality
  
  # From Relapse (only 1 cycle allowed) - goes to first tunnel state
  tpm["Relapse", "Subsequent remission_1"] <- 1 - mort_relapse
  tpm["Relapse", "Dead"]                   <- mort_relapse
  
  # From Subsequent remission tunnel states
  for (tunnel_idx in 1:max_tunnel_states) {
    tunnel_state <- paste0("Subsequent remission_", tunnel_idx)
    
    # Get mortality for this tunnel state
    mort_subseq <- get_steroid_ae_probs_by_state(
      cycle = cycle,
      state = "Subsequent remission",
      arm   = arm,
      duration = steroid_effect_duration_months,
      time_in_state = tunnel_idx  # tunnel_idx represents time in state
    )$mortality
    
    # Transitions from this tunnel state
    tpm[tunnel_state, "Relapse"] <- p_flare
    tpm[tunnel_state, "Dead"]    <- mort_subseq
    
    # Move to next tunnel state if not the last one
    if (tunnel_idx < max_tunnel_states) {
      next_tunnel <- paste0("Subsequent remission_", tunnel_idx + 1)
      tpm[tunnel_state, next_tunnel] <- 1 - p_flare - mort_subseq
    } else {
      # Stay in the last tunnel state (representing long-term subsequent remission)
      tpm[tunnel_state, tunnel_state] <- 1 - p_flare - mort_subseq
    }
  }
  
  # Absorbing state
  tpm["Dead", "Dead"] <- 1
  
  return(tpm)
}

# Initialize population arrays for tracking state populations
pop <- array(0, 
             dim = c(n_states, n_cycles + 1, n_arms),
             dimnames = list(state = state_names, 
                             cycle = 0:n_cycles, 
                             arm = arms))

# Initialize with everyone in Initial remission at cycle 0
pop["Initial remission", 1, ] <- 1

# Model GCA-related complications and steroid-related adverse events as events
event_log <- list(
  vl = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  stroke = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  mi = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  hf = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  fractures = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  diabetes = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  infection = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL))
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

# Sum costs and QALYs for each cycle
cycle_empty_array <- array(0,
                           dim = c(n_arms, n_cycles),
                           dimnames = list(arms = arms, cycle = NULL))

cycle_costs <- cycle_qalys <- cycle_empty_array
event_costs <- event_utility <- cycle_empty_array

# Running the model 
for (i in 1:n_arms) {
  for (j in 2:(n_cycles + 1)) {
    cycle_num <- j - 1  # Current cycle number (1-based)
    
    # Get TPM for this cycle
    tpm <- get_tpm(arm = i, cycle = cycle_num)
    
    # Update population using TPM
    pop[, j, i] <- pop[, j-1, i] %*% tpm
    
    # Process steroid-related adverse events
    # Define cycle time in years
    t_years <- cycle_num * cycle_length / 12
    current_year <- floor(t_years) + 1
    
    # Process Initial remission state
    if (pop["Initial remission", j, i] > 1e-10) {
      state_pop <- pop["Initial remission", j, i]
      
      steroid_ae_probs <- get_steroid_ae_probs_by_state(
        cycle = cycle_num,
        state = "Initial remission",
        arm = i,
        duration = steroid_effect_duration_months,
        time_in_state = cycle_num  # time_in_state = current cycle number
      )
      
      # Process steroid-related adverse events
      steroid_events <- c("stroke", "mi", "hf", "fractures", "diabetes", "infection")
      for (event in steroid_events) {
        if (event %in% names(steroid_ae_probs)) {
          prob_event <- steroid_ae_probs[[event]]
          n_event <- state_pop * prob_event
          event_log[[event]][i, cycle_num] <- event_log[[event]][i, cycle_num] + n_event
          
          # Update cumulative population
          if (cycle_num == 1) {
            cumulative_populations[[event]][i, cycle_num] <- n_event
          } else {
            cumulative_populations[[event]][i, cycle_num] <- cumulative_populations[[event]][i, cycle_num-1] + n_event
          }
        }
      }
    }
    
    # Process Relapse state
    if (pop_expanded["Relapse", j, i] > 1e-10) {
      state_pop <- pop_expanded["Relapse", j, i]
      
      # GCA-related complications: vision loss only occurs in relapse state for that cycle
      n_vl_incident <- state_pop * p_vl
      event_log$vl[i, cycle_num] <- event_log$vl[i, cycle_num] + n_vl_incident
      
      # Vision loss cases are only for this cycle (not cumulative)
      cumulative_populations$vl[i, cycle_num] <- n_vl_incident
      
      # Steroid-related adverse events in relapse
      steroid_ae_probs <- get_steroid_ae_probs_by_state(
        cycle = cycle_num,
        state = "Relapse",
        arm = i,
        duration = steroid_effect_duration_months,
        time_in_state = 1
      )
      
      for (event in steroid_events) {
        if (event %in% names(steroid_ae_probs)) {
          prob_event <- steroid_ae_probs[[event]]
          n_event <- state_pop * prob_event
          event_log[[event]][i, cycle_num] <- event_log[[event]][i, cycle_num] + n_event
          
          if (cycle_num == 1) {
            cumulative_populations[[event]][i, cycle_num] <- n_event
          } else {
            cumulative_populations[[event]][i, cycle_num] <- cumulative_populations[[event]][i, cycle_num-1] + n_event
          }
        }
      }
    }
    
    # Process Subsequent remission tunnel states
    for (tunnel_idx in 1:max_tunnel_states) {
      tunnel_state <- paste0("Subsequent remission_", tunnel_idx)
      
      if (pop[tunnel_state, j, i] > 1e-10) {
        state_pop <- pop[tunnel_state, j, i]
        
        # Steroid-related adverse events - tunnel_idx represents time_in_state
        steroid_ae_probs <- get_steroid_ae_probs_by_state(
          cycle = cycle_num,
          state = "Subsequent remission",
          arm = i,
          duration = steroid_effect_duration_months,
          time_in_state = tunnel_idx  # Each tunnel has different time_in_state
        )
        
        for (event in steroid_events) {
          if (event %in% names(steroid_ae_probs)) {
            prob_event <- steroid_ae_probs[[event]]
            n_event <- state_pop * prob_event
            event_log[[event]][i, cycle_num] <- event_log[[event]][i, cycle_num] + n_event
            
            if (cycle_num == 1) {
              cumulative_populations[[event]][i, cycle_num] <- n_event
            } else {
              cumulative_populations[[event]][i, cycle_num] <- cumulative_populations[[event]][i, cycle_num-1] + n_event
            }
          }
        }
      }
    }
    
    # Apply costs and utility impacts for all conditions
    # Vision loss - only for current cycle (no carry forward)
    if (cumulative_populations$vl[i, cycle_num] > 0) {
      vl_cost <- ifelse(current_year == 1, c_vl_y1, c_vl_y2on)
      event_costs[i, cycle_num] <- event_costs[i, cycle_num] + 
        apply_discount(cumulative_populations$vl[i, cycle_num] * vl_cost, disc_rate, t_years)
      event_utility[i, cycle_num] <- event_utility[i, cycle_num] + 
        apply_discount(cumulative_populations$vl[i, cycle_num] * utility_decrement_vl, disc_rate, t_years)
    }
    
    # Steroid-related adverse events costs and utility impacts
    for (event in steroid_events) {
      if (cumulative_populations[[event]][i, cycle_num] > 0) {
        # Event-specific cost and utility loss
        if (event == "stroke") {
          cost <- ifelse(current_year == 1, c_stroke_y1, c_stroke_y2on)
          qloss <- ifelse(current_year == 1, utility_decrement_stroke_y1, utility_decrement_stroke_y2on)
        } else if (event == "mi") {
          cost <- ifelse(current_year == 1, c_mi_y1, c_mi_y2on)
          qloss <- ifelse(current_year == 1, utility_decrement_mi_y1, utility_decrement_mi_y2on)
        } else if (event == "hf") {
          cost <- ifelse(current_year == 1, c_hf_y1, c_hf_y2on)
          qloss <- ifelse(current_year == 1, utility_decrement_hf_y1, utility_decrement_hf_y2on)
        } else if (event == "fractures") {
          cost <- ifelse(current_year == 1, c_fracture_y1, c_fracture_y2on)
          qloss <- ifelse(current_year == 1, utility_decrement_fracture_y1, utility_decrement_fracture_y2on)
        } else if (event == "diabetes") {
          cost <- c_dm
          qloss <- utility_decrement_dm
        } else if (event == "infection") {
          cost <- c_infection
          qloss <- q_infection
        } else {
          cost <- c_remission
          qloss <- q_remission
        }
        
        event_costs[i, cycle_num] <- event_costs[i, cycle_num] + 
          apply_discount(cumulative_populations[[event]][i, cycle_num] * cost, disc_rate, t_years)
        event_utility[i, cycle_num] <- event_utility[i, cycle_num] + 
          apply_discount(cumulative_populations[[event]][i, cycle_num] * qloss, disc_rate, t_years)
      }
    }
  }
}

# Create summary populations for the original 4 states for reporting
pop_summary <- array(0, 
                     dim = c(n_states, n_cycles + 1, n_arms),
                     dimnames = list(state = state_names, 
                                     cycle = 0:n_cycles, 
                                     arm = arms))

# Aggregate tunnel states back to main states
for (i in 1:n_arms) {
  for (j in 1:(n_cycles + 1)) {
    pop_summary["Initial remission", j, i] <- pop["Initial remission", j, i]
    pop_summary["Relapse", j, i] <- pop["Relapse", j, i]
    pop_summary["Subsequent remission", j, i] <- sum(pop[subseq_tunnel_names, j, i])
    pop_summary["Dead", j, i] <- pop["Dead", j, i]
  }
}