########################################
# Semi-Markov model with Tunnel States #
########################################

run_markov_model <- function(params = NULL, 
                             use_baseline = TRUE,
                             trans_probs_treat = trans_probs_first_flare_treat_exp,
                             trans_probs_control = trans_probs_first_flare_control_exp) {
  # params: named list of parameters to override baseline values
  # use_baseline: if TRUE, use baseline parameter values; if FALSE, must provide params
  # trans_probs_treat: transition probabilities for treatment arm (defaults to exponential)
  # trans_probs_control: transition probabilities for control arm (defaults to exponential)
  
  # If params provided, assign them to the environment
  if (!is.null(params)) {
    # Scalar parameters
    if ("scalar_params" %in% names(params)) {
      for (var_name in names(params$scalar_params)) {
        assign(var_name, params$scalar_params[[var_name]], envir = environment())
      }
    }
    
    # Vector parameters
    if ("vector_params" %in% names(params)) {
      for (vec_name in names(params$vector_params)) {
        assign(vec_name, params$vector_params[[vec_name]], envir = environment())
      }
    }
  }
  
  # If use_baseline = TRUE and params not provided, the function will use
  # whatever parameter values are in the global environment

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
    # Check if we're using sampled probabilities (PSA) or baseline
    if (is.data.frame(trans_probs_treat)) {
      # Baseline mode: trans_probs is a data frame
      p_first_flare <- trans_probs_treat$trans_prob[cycle]
    } else {
      # PSA mode: trans_probs is a vector of probabilities for this simulation
      p_first_flare <- trans_probs_treat[cycle]
    }
    p_flare <- probs_subseq_flare_treat
  } else {
    # Check if we're using sampled probabilities (PSA) or baseline
    if (is.data.frame(trans_probs_control)) {
      # Baseline mode: trans_probs is a data frame
      p_first_flare <- trans_probs_control$trans_prob[cycle]
    } else {
      # PSA mode: trans_probs is a vector of probabilities for this simulation
      p_first_flare <- trans_probs_control[cycle]
    }
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

# Track new incidence of GCA-related complications and steroid-related adverse events 
event_log <- list(
  vl = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  stroke = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  mi = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  hf = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  fractures = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  dm = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL)),
  infection = array(0, dim = c(n_arms, n_cycles), dimnames = list(arms = arms, cycle = NULL))
)

# Track total number of patients who have ever experienced a specific event up to a given cycle
cumulative_populations <- list(
  vl = array(0, dim = c(n_arms, n_cycles + 1), dimnames = list(arms = arms, cycle = NULL)),
  stroke = array(0, dim = c(n_arms, n_cycles + 1), dimnames = list(arms = arms, cycle = NULL)),
  mi = array(0, dim = c(n_arms, n_cycles + 1), dimnames = list(arms = arms, cycle = NULL)),
  hf = array(0, dim = c(n_arms, n_cycles + 1), dimnames = list(arms = arms, cycle = NULL)),
  fractures = array(0, dim = c(n_arms, n_cycles + 1), dimnames = list(arms = arms, cycle = NULL)),
  dm = array(0, dim = c(n_arms, n_cycles + 1), dimnames = list(arms = arms, cycle = NULL)),
  infection = array(0, dim = c(n_arms, n_cycles + 1), dimnames = list(arms = arms, cycle = NULL))
)

# Sum costs and QALYs for each cycle
cycle_empty_array <- array(0,
                           dim = c(n_arms, n_cycles),
                           dimnames = list(arms = arms, cycle = NULL))

cycle_costs <- cycle_qalys <- cycle_empty_array
event_costs <- event_qalys <- cycle_empty_array

# Track costs and QALY impacts separately by event 
all_events <- c("vl", "stroke", "mi", "hf", "fractures", "dm", "infection")

event_costs_by_event <- array(
  0,
  dim = c(length(all_events), n_arms, n_cycles),
  dimnames = list(event = all_events, arm = arms, cycle = NULL)
)

event_qalys_by_event <- array(
  0,
  dim = c(length(all_events), n_arms, n_cycles),
  dimnames = list(event = all_events, arm = arms, cycle = NULL)
)

# Function to determine if patient is on steroid treatment
is_on_steroid <- function(cycle, state, arm, time_in_state = 1) {
  if (state == "Initial remission") {
    if (arm == 1) {  # Treatment group (Tocilizumab)
      return(cycle <= 6)
    } else {  # Control group
      return(cycle <= 14)
    }
  } else if (state == "Relapse") {
    return(TRUE)  # Always on steroid during relapse (only 1 cycle)
  } else if (state == "Subsequent remission") {
    if (arm == 1) {  # Treatment group (Tocilizumab)
      return(time_in_state <= 2)
    } else {  # Control group
      return(time_in_state <= 11)
    }
  }
  return(FALSE)
}

# Function to determine if this is the start of a steroid treatment episode
is_steroid_start <- function(cycle, state, arm, time_in_state = 1) {
  if (state == "Initial remission") {
    # Start of initial steroid treatment
    return(cycle == 1)
  } else if (state == "Relapse") {
    # Always start steroid treatment when entering relapse
    return(TRUE)  # Since relapse is only 1 cycle, this is always the start
  } else if (state == "Subsequent remission") {
    # Start of subsequent steroid treatment (first tunnel state)
    return(time_in_state == 1)
  }
  return(FALSE)
}

# Running the model 
for (i in 1:n_arms) {
  for (j in 2:(n_cycles + 1)) {
    cycle_num <- j - 1  # Current cycle number (1-based)
    
    # Get TPM for this cycle
    tpm <- get_tpm(arm = i, cycle = cycle_num)
    
    # Update population using TPM
    pop[, j, i] <- pop[, j-1, i] %*% tpm
    
    # Calculate costs and utilities for this cycle
    t_years <- cycle_num * cycle_length / 12
    current_year <- floor(t_years) + 1
    
    # Reset cycle-specific event log before calculations
    for (event in c("vl", steroid_events)) {
      event_log[[event]][i, cycle_num] <- 0
    }
    
    # Process Initial remission state
    if (pop["Initial remission", j, i] > 1e-10) {
      state_pop <- pop["Initial remission", j, i]
      on_steroid <- is_on_steroid(cycle_num, "Initial remission", i, cycle_num)
      steroid_start <- is_steroid_start(cycle_num, "Initial remission", i, cycle_num)
      
      # State costs
      state_cost <- 0
      state_utility <- q_remission
      
      # Treatment costs
      if (i == 1) {  # Treatment arm
        if (is_on_tcz(cycle_num)) { # TCZ monthly cost during first 12 cycles only
          state_cost <- state_cost + c_TCZ + c_TCZ_monitoring
        }
        if (on_steroid) {
          state_utility <- state_utility - q_decrement_steroid
          if (steroid_start) {
            state_cost <- state_cost + c_steroid_treat_init_remission 
          }
        }
      } else {  # Control arm
        if (cycle_num<=12) { # steroid monitoring cost assumed to occur during first 12 cycles only as well
          state_cost <- state_cost + c_steroid_monitoring
        }
        if (on_steroid) {
          state_utility <- state_utility - q_decrement_steroid
          if (steroid_start) {
            state_cost <- state_cost + c_steroid_control_init_remission 
          }
        }
      }
      
      # Apply discounting and add to cycle totals
      cycle_costs[i, cycle_num] <- cycle_costs[i, cycle_num] + 
        apply_discount(state_pop * state_cost, disc_rate, t_years)
      cycle_qalys[i, cycle_num] <- cycle_qalys[i, cycle_num] + 
        apply_discount(state_pop * state_utility, disc_rate, t_years)
      
      # Process steroid-related adverse events
      steroid_ae_probs <- get_steroid_ae_probs_by_state(
        cycle = cycle_num,
        state = "Initial remission",
        arm = i,
        duration = steroid_effect_duration_months,
        time_in_state = cycle_num
      )
      
      for (event in steroid_events) {
        if (event %in% names(steroid_ae_probs)) {
          prob_event <- steroid_ae_probs[[event]]
          
          # Calculate new incident cases 
          n_event <- pop["Initial remission", j-1, i] * prob_event
          
          # Update incident event log: accumulating incident events in multiple states
          event_log[[event]][i, cycle_num] <- event_log[[event]][i, cycle_num] + n_event
          
          }
        }
      }
    
    # Process Relapse state
    if (pop["Relapse", j, i] > 1e-10) {
      state_pop <- pop["Relapse", j, i]
      on_steroid <- is_on_steroid(cycle_num, "Relapse", i, 1)
      steroid_start <- is_steroid_start(cycle_num, "Relapse", i, 1)
      
      # State costs and utilities
      state_cost <- 0
      state_utility <- q_relapse
      
      # Treatment costs
      if (i == 1) {  # Treatment arm
        if (is_on_tcz(cycle_num)) { # TCZ monthly cost during first 12 cycles only
          state_cost <- state_cost + c_TCZ + c_TCZ_monitoring
        }
        if (on_steroid) {
          state_utility <- state_utility - q_decrement_steroid
          if (steroid_start) {
            state_cost <- state_cost + c_steroid_treat_relapse 
          }
        }
      } else {  # Control arm
        if (cycle_num<=12) { # steroid monitoring cost assumed to occur during first 12 cycles only as well
          state_cost <- state_cost + c_steroid_monitoring
        }
        if (on_steroid) {
          state_utility <- state_utility - q_decrement_steroid
          if (steroid_start) {
            state_cost <- state_cost + c_steroid_treat_relapse 
          }
        }
      }
      
      # Apply discounting and add to cycle totals
      # Note: The population for costs should be the population in the state *during* the cycle (pop[, j, i])
      cycle_costs[i, cycle_num] <- cycle_costs[i, cycle_num] +
        apply_discount(state_pop * state_cost, disc_rate, t_years)
      cycle_qalys[i, cycle_num] <- cycle_qalys[i, cycle_num] +
        apply_discount(state_pop * state_utility, disc_rate, t_years)
      
      # GCA-related complications (Vision Loss)
      if (pop["Relapse", j-1, i] > 1e-10) {
        n_vl_incident <- pop["Relapse", j-1, i] * p_vl
        # Update the incident event log by adding to the existing count for this cycle
        event_log$vl[i, cycle_num] <- event_log$vl[i, cycle_num] + n_vl_incident
      }
      
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
          
          # Calculate new incident cases 
          n_event <- pop["Relapse", j-1, i] * prob_event
          
          # Update incident event log: accumulating incident events in multiple states
          event_log[[event]][i, cycle_num] <- event_log[[event]][i, cycle_num] + n_event
        }
      }
    }
    
    # Process Subsequent remission tunnel states
    for (tunnel_idx in 1:max_tunnel_states) {
      tunnel_state <- paste0("Subsequent remission_", tunnel_idx)
      
      if (pop[tunnel_state, j, i] > 1e-10) {
        state_pop <- pop[tunnel_state, j, i]
        on_steroid <- is_on_steroid(cycle_num, "Subsequent remission", i, tunnel_idx)
        steroid_start <- is_steroid_start(cycle_num, "Subsequent remission", i, tunnel_idx)
        
        # State costs and utilities
        state_cost <- 0
        state_utility <- q_remission
        
        # Treatment costs
        if (i == 1) {  # Treatment arm
          if (is_on_tcz(cycle_num)) { # TCZ monthly cost during first 12 cycles only
            state_cost <- state_cost + c_TCZ + c_TCZ_monitoring
          }
          if (on_steroid) {
            state_utility <- state_utility - q_decrement_steroid
            if (steroid_start) {
              state_cost <- state_cost + c_steroid_treat_sub_remission 
            }
          }
        } else {  # Control arm
          if (cycle_num<=12) { # steroid monitoring cost assumed to occur during first 12 cycles only as well
            state_cost <- state_cost + c_steroid_monitoring
          }
          if (on_steroid) {
            state_utility <- state_utility - q_decrement_steroid
            if (steroid_start) {
              state_cost <- state_cost + c_steroid_control_sub_remission 
            }
          }
        }
        
        # Apply discounting and add to cycle totals
        # Use state_pop_for_costs here
        cycle_costs[i, cycle_num] <- cycle_costs[i, cycle_num] + 
          apply_discount(state_pop * state_cost, disc_rate, t_years)
        cycle_qalys[i, cycle_num] <- cycle_qalys[i, cycle_num] + 
          apply_discount(state_pop * state_utility, disc_rate, t_years)
        
        # Process steroid-related adverse events
        steroid_ae_probs <- get_steroid_ae_probs_by_state(
          cycle = cycle_num,
          state = "Subsequent remission",
          arm = i,
          duration = steroid_effect_duration_months,
          time_in_state = tunnel_idx
        )
        
        for (event in steroid_events) {
          if (event %in% names(steroid_ae_probs)) {
            prob_event <- steroid_ae_probs[[event]]
            
            # Calculate new incident cases from the population at risk
            if (pop[tunnel_state, j-1, i] > 1e-10) {
              n_event <- pop[tunnel_state, j-1, i] * prob_event
              
              # Update incident event log
              event_log[[event]][i, cycle_num] <- event_log[[event]][i, cycle_num] + n_event
            }
          }
        }
      }
    }
    
    # Update the cumulative_populations arrays once per cycle, after all states have been processed.
    for (event in all_events) {
      # The cumulative population is the total from the last cycle plus the new incidents
      if (cycle_num == 1) {
        cumulative_populations[[event]][i, cycle_num] <- event_log[[event]][i, cycle_num]
      } else {
        cumulative_populations[[event]][i, cycle_num] <- cumulative_populations[[event]][i, cycle_num-1] + event_log[[event]][i, cycle_num]
      }
    }
    
    # Apply costs and utility impacts for all adverse events
    for (event in all_events) {
      
      # YEAR 1: incident cases this cycle 
      if (event_log[[event]][i, cycle_num] > 1e-10) {
        new_cases <- event_log[[event]][i, cycle_num]
        
        # Get year 1 costs and QALY decrements for new cases
        if (event == "vl") {
          cost_y1 <- c_vl_y1/12
          decrement_y1 <- (q_remission - q_vl)/12
        } else if (event == "stroke") {
          cost_y1 <- c_stroke_y1/12
          decrement_y1 <- (q_remission - q_stroke_y1)/12
        } else if (event == "mi") {
          cost_y1 <- c_mi_y1/12
          decrement_y1 <- (q_remission - q_mi_y1)/12      
        } else if (event == "hf") {
          cost_y1 <- c_hf_y1/12
          decrement_y1 <- (q_remission - q_hf_y1)/12      
        } else if (event == "fractures") {
          cost_y1 <- c_fracture_y1/12
          decrement_y1 <- (q_remission - q_fracture_y1)/12 
        } else if (event == "dm") {
          cost_y1 <- c_dm/12
          decrement_y1 <- (q_remission - q_dm)/12        
        } else if (event == "infection") {
          cost_y1 <- c_infection/12
          decrement_y1 <- q_decrement_infection/12  
        } else {
          cost_y1 <- 0
          decrement_y1 <- 0
        }
     
        # Apply year 1 costs and utilities to new cases
        event_costs[i, cycle_num] <- event_costs[i, cycle_num] + 
          apply_discount(new_cases * cost_y1, disc_rate, t_years)
        event_qalys[i, cycle_num] <- event_qalys[i, cycle_num] - 
          apply_discount(new_cases * decrement_y1, disc_rate, t_years)
        
        # Per-event tracking
        idx <- match(event, all_events)
        event_costs_by_event[idx, i, cycle_num] <- event_costs_by_event[idx, i, cycle_num] +
          apply_discount(new_cases * cost_y1, disc_rate, t_years)
        event_qalys_by_event[idx, i, cycle_num] <- event_qalys_by_event[idx, i, cycle_num] -
          apply_discount(new_cases * decrement_y1, disc_rate, t_years)
      }
      
      # YEARS 2 ONWARDS: Apply year 2+ costs/utilities to all existing cases from previous cycles
      if (cycle_num > 1) {
        # Get the cumulative population up to the end of the previous cycle
        existing_cumulative <- cumulative_populations[[event]][i, j-1]
        
        if (existing_cumulative > 1e-10) {
          # Get year 2+ costs and QALY decrements for existing cases
          if (event == "vl") {
            cost_y2on <- c_vl_y2on/12
            decrement_y2on <- (q_remission - q_vl)/12
          } else if (event == "stroke") {
            cost_y2on <- c_stroke_y2on/12
            decrement_y2on <- (q_remission - q_stroke_y2on)/12  
          } else if (event == "mi") {
            cost_y2on <- c_mi_y2on/12
            decrement_y2on <- (q_remission - q_mi_y2on)/12      
          } else if (event == "hf") {
            cost_y2on <- c_hf_y2on/12
            decrement_y2on <- (q_remission - q_hf_y2on)/12      
          } else if (event == "fractures") {
            cost_y2on <- c_fracture_y2on/12
            decrement_y2on <- (q_remission - q_fracture_y2on)/12
          } else if (event == "dm") {
            cost_y2on <- c_dm/12
            decrement_y2on <- (q_remission - q_dm)/12            
          } else if (event == "infection") {
            cost_y2on <- 0 # No ongoing cost for infections
            decrement_y2on <- 0  
          } else {
            cost_y2on <- 0
            decrement_y2on <- 0
          }
          
          # Apply year 2+ costs and utilities to existing cases
          event_costs[i, cycle_num] <- event_costs[i, cycle_num] + 
            apply_discount(existing_cumulative * cost_y2on, disc_rate, t_years)
          event_qalys[i, cycle_num] <- event_qalys[i, cycle_num] - 
            apply_discount(existing_cumulative * decrement_y2on, disc_rate, t_years)
          
          # Per-event tracking
          idx <- matchidx <- match(event, all_events)
          event_costs_by_event[idx, i, cycle_num] <- event_costs_by_event[idx, i, cycle_num] +
            apply_discount(existing_cumulative * cost_y2on, disc_rate, t_years)
          event_qalys_by_event[idx, i, cycle_num] <- event_qalys_by_event[idx, i, cycle_num] -
            apply_discount(existing_cumulative * decrement_y2on, disc_rate, t_years)
        }
      }
    }
    
  }
}

# Calculate total costs and QALYs
total_costs <- cycle_costs + event_costs
total_qalys <- cycle_qalys + event_qalys

# Create summary populations for the original 4 states for reporting (need to update)
pop_summary <- array(0, 
                     dim = c(4, n_cycles + 1, n_arms),
                     dimnames = list(state = c("Initial remission", "Relapse", "Subsequent remission", "Dead"), 
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

# Summarize total outcomes
total_costs_sum <- rowSums(total_costs)
total_qalys_sum <- rowSums(total_qalys)
inc_cost  <- total_costs_sum[1] - total_costs_sum[2]
inc_qalys <- total_qalys_sum[1] - total_qalys_sum[2]
icer <- inc_cost / inc_qalys

# Return all results as a list
return(list(
  total_costs_sum = total_costs_sum,
  total_qalys_sum = total_qalys_sum,
  inc_cost = inc_cost,
  inc_qalys = inc_qalys,
  icer = icer,
  pop = pop,
  pop_summary = pop_summary,
  cycle_costs = cycle_costs,
  cycle_qalys = cycle_qalys,
  event_costs = event_costs,
  event_qalys = event_qalys,
  event_log = event_log,
  cumulative_populations = cumulative_populations,
  event_costs_by_event = event_costs_by_event,
  event_qalys_by_event = event_qalys_by_event
))
}

# SUMMARY OF ADVERSE EVENTS BY TREATMENT GROUP
################################################

summarize_adverse_events <- function(results) {
  
  # Extract necessary data from results
  event_log <- results$event_log
  pop_summary <- results$pop_summary
  event_costs_by_event <- results$event_costs_by_event
  event_qalys_by_event <- results$event_qalys_by_event

# Calculate total number of INCIDENT events over all cycles (for costing purposes)
total_incident_events_summary <- data.frame(
  Event = character(),
  TCZ_plus_GC = numeric(),
  GC_alone = numeric(),
  Difference = numeric(),
  stringsAsFactors = FALSE
)

# Vision loss events (sum of incident events)
total_vl_tcz <- sum(event_log$vl[1, ])
total_vl_control <- sum(event_log$vl[2, ])
total_incident_events_summary <- rbind(total_incident_events_summary, 
                                       data.frame(Event = "Vision Loss", 
                                                  TCZ_plus_GC = round(total_vl_tcz, 4),
                                                  GC_alone = round(total_vl_control, 4),
                                                  Difference = round(total_vl_tcz - total_vl_control, 4)))

# Steroid-related adverse events (sum of incident events)
steroid_events <- c("stroke", "mi", "hf", "fractures", "dm", "infection")
event_labels <- c("Stroke", "Myocardial Infarction", "Heart Failure", "Fractures", "Diabetes", "Infection")

for (i in seq_along(steroid_events)) {
  event <- steroid_events[i]
  label <- event_labels[i]
  
  total_event_tcz <- sum(event_log[[event]][1, ])
  total_event_control <- sum(event_log[[event]][2, ])
  
  total_incident_events_summary <- rbind(total_incident_events_summary, 
                                         data.frame(Event = label,
                                                    TCZ_plus_GC = round(total_event_tcz, 3),
                                                    GC_alone = round(total_event_control, 3),
                                                    Difference = round(total_event_tcz - total_event_control, 3)))
}


# Calculate event rates per 100 patient-years for each arm
alive_pop <- array(0, dim = c(n_arms, n_cycles))
for (i in 1:n_arms) {
  for (j in 1:n_cycles) {
    # Population alive at the start of cycle j
    alive_pop[i, j] <- 1 - pop_summary["Dead", j, i]  # Using j instead of j+1 since we want start of cycle
  }
}
total_patient_years <- rowSums(alive_pop) / 12  # Convert months to years

event_rates_per_100py_summary <- total_incident_events_summary 
event_rates_per_100py_summary$TCZ_plus_GC <- round((event_rates_per_100py_summary$TCZ_plus_GC / total_patient_years[1]) * 100, 2)
event_rates_per_100py_summary$GC_alone <- round((event_rates_per_100py_summary$GC_alone / total_patient_years[2]) * 100, 2)
event_rates_per_100py_summary$Difference <- event_rates_per_100py_summary$TCZ_plus_GC - event_rates_per_100py_summary$GC_alone

# Lifetime per-event cost and QALY impact summaries
# Sum over cycles -> [event x arm] totals
total_event_costs_by_event <- apply(event_costs_by_event, c(1, 2), sum)
total_event_qalys_by_event <- apply(event_qalys_by_event, c(1, 2), sum)  # negative values = QALY loss

# Build a tidy summary table
per_event_summary <- data.frame(
  Event      = total_incident_events_summary[,"Event"],
  TCZ_plus_GC_cost = round(total_event_costs_by_event[, "1"], 2),
  GC_alone_cost    = round(total_event_costs_by_event[, "0"], 2),
  Difference_cost  = round(total_event_costs_by_event[, "1"] - total_event_costs_by_event[, "0"], 2),
  TCZ_plus_GC_QALY = round(total_event_qalys_by_event[, "1"], 4),
  GC_alone_QALY    = round(total_event_qalys_by_event[, "0"], 4),
  Difference_QALY  = round(total_event_qalys_by_event[, "1"] - total_event_qalys_by_event[, "0"], 4),
  check.names = FALSE
)

return(list(
  total_incident_events = total_incident_events_summary,
  event_rates_per_100py = event_rates_per_100py_summary,
  total_patient_years = total_patient_years,
  per_event_cost_qaly = per_event_summary
))

}
