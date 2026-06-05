########################################
# Semi-Markov model with Tunnel States #
########################################

run_markov_model <- function(params = NULL, 
                             use_baseline = TRUE,
                             subgroup = "itt",
                             initial_age = 71,
                             time_horizon_months = 30*12,
                             time_horizon_weeks = 30*52, 
                             treatment_duration_months = 12,
                             steroid_effect_duration_months = 30*12) {
  
  # params: named list of parameters to override baseline values
  # use_baseline: if TRUE, use baseline parameter values; if FALSE, must provide params
  # subgroup: one of "newlydiag", "relapsed", or "itt" to determine which subsequent flare probabilities to use
  # initial_age: starting age of cohort
  # time_horizon_months: total time horizon in months
  # time_horizon_weeks: total time horizon in weeks
  # treatment_duration_months: duration of treatment in months
  # steroid_effect_duration_months: duration in months that steroid adverse event risks apply after starting steroids
  
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
  
  # If use_baseline = TRUE and params not provided, the function will use whatever parameter values are in the global environment
  
  # For "Subsequent remission", we need tunnel states to track time-in-state
  # Create tunnel states for subsequent remission (up to maximum duration needed)
  max_tunnel_states <- steroid_effect_duration_months  # Track up to steroid effect duration
  subseq_tunnel_names <- paste0("Subsequent remission_", 1:max_tunnel_states)
  
  # Number of states
  state_names <- c("Initial remission", "Relapse", subseq_tunnel_names, "Dead")
  n_states <- length(state_names)
  
  # Get patients' age and mortality rate in each cycle
  v.age <- initial_age + ((1: time_horizon_months - 1) / 12) 
  v.allcause.mort.monthly <- get_monthly_mort(v.age, v.allcause.mort)
  
  # Convert annual probabilities into monthly probabilities
  annual_probs <- c("p_vl", "p_stroke", "p_mi", "p_hf", "p_fractures", "p_infection", "p_dm") 
  
  for (prob in annual_probs) {
    monthly_prob <- 1 - exp(-(-log(1 - get(prob))) / 12)
    assign(paste0(prob, "_monthly"), monthly_prob)
  }
  
  # Function to calculate probabilities for steroid adverse events
  get_steroid_ae_probs_by_state <- function(cycle, state, arm, 
                                            time_in_state = 1,
                                            duration = steroid_effect_duration_months) {
    
    # Choose which steroid tapering schedule to use
    if (state == "Initial remission") {
      schedule_index <- ifelse(arm == 1, "treat", "control")  
    } else if (state %in% c("Relapse","Subsequent remission")) {
      schedule_index <- ifelse(arm == 1, "treat.relapse", "control.relapse")
    }
    
    max_cycle <- nrow(steroid.hr)
    
    # Determine index into HR vectors
    if (cycle > duration) {
      # steroid effect has worn off
      steroid_ae_probs <- list(
        stroke     = p_stroke_monthly,
        mi         = p_mi_monthly,
        hf         = p_hf_monthly,
        dm         = p_dm_monthly,
        infection  = p_infection_monthly,
        fractures  = p_fractures_monthly,
        mortality  = v.allcause.mort.monthly[cycle]
      )
      
    } else {
      if (state == "Initial remission") {
        idx <- min(cycle, max_cycle)
      } else if (state == "Relapse") {
        idx <- 1
      } else if (state == "Subsequent remission") {
        idx <- min(time_in_state + 1, max_cycle)
      }
      
      # Select HRs by schedule
      if (schedule_index == "treat") {
        cvd_hr <- v.steroid.hr.cvd.treat[idx]
        dm_hr <- v.steroid.hr.diabetes.treat[idx]
        infection_hr <- v.steroid.hr.infection.treat[idx]
        fractures_hr <- v.steroid.hr.fractures.treat[idx]
        mortality_hr <- v.steroid.hr.mortality.treat[idx]
      } else if (schedule_index == "control") {
        cvd_hr <- v.steroid.hr.cvd.control[idx]
        dm_hr <- v.steroid.hr.diabetes.control[idx]
        infection_hr <- v.steroid.hr.infection.control[idx]
        fractures_hr <- v.steroid.hr.fractures.control[idx]
        mortality_hr <- v.steroid.hr.mortality.control[idx]
      } else if (schedule_index == "treat.relapse") {
        cvd_hr <- v.steroid.hr.cvd.treat.relapse[idx]
        dm_hr <- v.steroid.hr.diabetes.treat.relapse[idx]
        infection_hr <- v.steroid.hr.infection.treat.relapse[idx]
        fractures_hr <- v.steroid.hr.fractures.treat.relapse[idx]
        mortality_hr <- v.steroid.hr.mortality.treat.relapse[idx]
      } else if (schedule_index == "control.relapse") {
        cvd_hr <- v.steroid.hr.cvd.control.relapse[idx]
        dm_hr <- v.steroid.hr.diabetes.control.relapse[idx]
        infection_hr <- v.steroid.hr.infection.control.relapse[idx]
        fractures_hr <- v.steroid.hr.fractures.control.relapse[idx]
        mortality_hr <- v.steroid.hr.mortality.control.relapse[idx]
      }
      
      steroid_ae_probs <- list(
        stroke     = p_stroke_monthly    * cvd_hr,
        mi         = p_mi_monthly        * cvd_hr,
        hf         = p_hf_monthly        * cvd_hr,
        dm         = p_dm_monthly        * dm_hr,
        infection  = p_infection_monthly * infection_hr,
        fractures  = p_fractures_monthly * fractures_hr,
        mortality  = v.allcause.mort.monthly[cycle] * mortality_hr
      )
    }
    
    return(steroid_ae_probs)
  }
  
  # Time-dependent transition probability matrices (TPM) creation function
  get_tpm <- function(arm, cycle) {
    
    # Ensure cycle does not exceed the length of input data
    cycle <- min(cycle, time_horizon_months) 
    
    # Extract cycle-specific flare probabilities
    if (arm == 1) {
      # Check if we're using sampled probabilities (PSA) or baseline
      if (is.data.frame(trans_probs_first_flare_treat)) {
        # Baseline mode: trans_probs is a data frame
        p_first_flare <- trans_probs_first_flare_treat$trans_prob[cycle]
      } else {
        # PSA mode: trans_probs is a vector of probabilities for this simulation
        p_first_flare <- trans_probs_first_flare_treat[cycle]
      }
      
      # Select appropriate subsequent flare probability based on subgroup
      if (subgroup == "newlydiag") {
        p_flare <- probs_subseq_flare_newlydiag_treat
      } else if (subgroup == "relapsed") {
        p_flare <- probs_subseq_flare_relapsed_treat
      } else if (subgroup == "itt") {
        p_flare <- probs_subseq_flare_itt_treat
      } else {
        stop("Invalid subgroup. Must be one of: 'newlydiag', 'relapsed', or 'itt'")
      }
      
    } else {
      # Check if we're using sampled probabilities (PSA) or baseline
      if (is.data.frame(trans_probs_first_flare_control)) {
        # Baseline mode: trans_probs is a data frame
        p_first_flare <- trans_probs_first_flare_control$trans_prob[cycle]
      } else {
        # PSA mode: trans_probs is a vector of probabilities for this simulation
        p_first_flare <- trans_probs_first_flare_control[cycle]
      }
      
      # Select appropriate subsequent flare probability based on subgroup
      if (subgroup == "newlydiag") {
        p_flare <- probs_subseq_flare_newlydiag_control
      } else if (subgroup == "relapsed") {
        p_flare <- probs_subseq_flare_relapsed_control
      } else if (subgroup == "itt") {
        p_flare <- probs_subseq_flare_itt_control
      } else {
        stop("Invalid subgroup. Must be one of: 'newlydiag', 'relapsed', or 'itt'")
      }
      
    }
    
    # Initialize TPM
    tpm <- matrix(0, nrow = n_states, ncol = n_states,
                  dimnames = list(from = state_names, to = state_names))
    
    # Get mortality for Initial remission 
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
  
  # Apply costs and utility impacts for all adverse events
  get_event_params <- function(event) {
    params <- switch(event,
                     "vl" = list(
                       cost_y1 = c_vl_y1 / 12,
                       cost_y2on = c_vl_y2on / 12,
                       decrement_y1 = (q_remission - q_vl) / 12,
                       decrement_y2on = (q_remission - q_vl) / 12
                     ),
                     "stroke" = list(
                       cost_y1 = c_stroke_y1 / 12,
                       cost_y2on = c_stroke_y2on / 12,
                       decrement_y1 = (q_remission - q_stroke_y1) / 12,
                       decrement_y2on = (q_remission - q_stroke_y2on) / 12
                     ),
                     "mi" = list(
                       cost_y1 = c_mi_y1 / 12,
                       cost_y2on = c_mi_y2on / 12,
                       decrement_y1 = (q_remission - q_mi_y1) / 12,
                       decrement_y2on = (q_remission - q_mi_y2on) / 12
                     ),
                     "hf" = list(
                       cost_y1 = c_hf_y1 / 12,
                       cost_y2on = c_hf_y2on / 12,
                       decrement_y1 = (q_remission - q_hf_y1) / 12,
                       decrement_y2on = (q_remission - q_hf_y2on) / 12
                     ),
                     "fractures" = list(
                       cost_y1 = c_fracture_y1 / 12,
                       cost_y2on = c_fracture_y2on / 12,
                       decrement_y1 = (q_remission - q_fracture_y1) / 12,
                       decrement_y2on = (q_remission - q_fracture_y2on) / 12
                     ),
                     "dm" = list(
                       cost_y1 = c_dm / 12,
                       cost_y2on = c_dm / 12,
                       decrement_y1 = (q_remission - q_dm) / 12,
                       decrement_y2on = (q_remission - q_dm) / 12
                     ),
                     "infection" = list(
                       cost_y1 = c_infection / 12,
                       cost_y2on = 0,  # No ongoing cost
                       decrement_y1 = q_decrement_infection / 12,
                       decrement_y2on = 0  # No ongoing utility decrement
                     ),
                     list(cost_y1 = 0, cost_y2on = 0, decrement_y1 = 0, decrement_y2on = 0)
    )
    return(params)
  }
  
  # Initialize population arrays for tracking state populations
  pop <- array(0, 
               dim = c(n_states, time_horizon_months + 1, n_arms),
               dimnames = list(state = state_names, 
                               cycle = 0:time_horizon_months, 
                               arm = arms))
  
  # Initialize with everyone in Initial remission at cycle 0
  pop["Initial remission", 1, ] <- 1
  
  # Track new incidence of GCA-related complications and steroid-related adverse events 
  event_log <- list(
    vl = array(0, dim = c(n_arms, time_horizon_months), dimnames = list(arms = arms, cycle = NULL)),
    stroke = array(0, dim = c(n_arms, time_horizon_months), dimnames = list(arms = arms, cycle = NULL)),
    mi = array(0, dim = c(n_arms, time_horizon_months), dimnames = list(arms = arms, cycle = NULL)),
    hf = array(0, dim = c(n_arms, time_horizon_months), dimnames = list(arms = arms, cycle = NULL)),
    fractures = array(0, dim = c(n_arms, time_horizon_months), dimnames = list(arms = arms, cycle = NULL)),
    dm = array(0, dim = c(n_arms, time_horizon_months), dimnames = list(arms = arms, cycle = NULL)),
    infection = array(0, dim = c(n_arms, time_horizon_months), dimnames = list(arms = arms, cycle = NULL))
  )
  
  # Track total number of patients who have ever experienced a specific event up to a given cycle
  cumulative_populations <- list(
    vl = array(0, dim = c(n_arms, time_horizon_months + 1), dimnames = list(arms = arms, cycle = NULL)),
    stroke = array(0, dim = c(n_arms, time_horizon_months + 1), dimnames = list(arms = arms, cycle = NULL)),
    mi = array(0, dim = c(n_arms, time_horizon_months + 1), dimnames = list(arms = arms, cycle = NULL)),
    hf = array(0, dim = c(n_arms, time_horizon_months + 1), dimnames = list(arms = arms, cycle = NULL)),
    fractures = array(0, dim = c(n_arms, time_horizon_months + 1), dimnames = list(arms = arms, cycle = NULL)),
    dm = array(0, dim = c(n_arms, time_horizon_months + 1), dimnames = list(arms = arms, cycle = NULL)),
    infection = array(0, dim = c(n_arms, time_horizon_months + 1), dimnames = list(arms = arms, cycle = NULL))
  )
  
  # Sum costs and QALYs for each cycle
  cycle_empty_array <- array(0,
                             dim = c(n_arms, time_horizon_months),
                             dimnames = list(arms = arms, cycle = NULL))
  
  cycle_costs <- cycle_qalys <- cycle_empty_array
  event_costs <- event_qalys <- cycle_empty_array
  
  # Track costs and QALY impacts separately by event 
  all_events <- c("vl", "stroke", "mi", "hf", "fractures", "dm", "infection")
  
  event_costs_by_event <- array(
    0,
    dim = c(length(all_events), n_arms, time_horizon_months),
    dimnames = list(event = all_events, arm = arms, cycle = NULL)
  )
  
  event_qalys_by_event <- array(
    0,
    dim = c(length(all_events), n_arms, time_horizon_months),
    dimnames = list(event = all_events, arm = arms, cycle = NULL)
  )
  
  
  # Running the model 
  for (i in 1:n_arms) {
    for (j in 2:(time_horizon_months + 1)) {
      cycle_num <- j - 1  # Current cycle number (1-based)
      
      # Get TPM for this cycle
      tpm <- get_tpm(arm = i, cycle = cycle_num)
      
      # Update population using TPM
      pop[, j, i] <- pop[, j-1, i] %*% tpm
      
      # Calculate costs and utilities for this cycle
      t_years <- (cycle_num - 1 ) / 12
      current_year <- floor(t_years) + 1
      
      # Process Initial remission state
      if (pop["Initial remission", j, i] > 1e-10) {
        state_pop <- pop["Initial remission", j, i]
        
        # State costs and utilities
         if (i == 1) {  # Treatment arm
           
          if (cycle_num <= treatment_duration_months) { # TCZ monthly cost during treatment duration only
            state_cost <- ifelse(cycle_num == 1, c_steroid_treat_init_remission, 0) + c_TCZ + c_TCZ_monitoring
          } else {
            state_cost <- 0
          }
           
          if (cycle_num <= 6) {  # On steroids for 6 months in treatment arm
            state_utility <- (q_remission - q_decrement_steroid)/12
          } else {
            state_utility <- q_remission/12
          }
           
        } else {  # Control arm
          
          if (cycle_num <= 13) {  # On steroids for about 13 months 
            state_cost <- ifelse(cycle_num == 1, c_steroid_control_init_remission, 0) + c_steroid_monitoring
            state_utility <- (q_remission - q_decrement_steroid)/12
          } else {
            state_cost <- 0
            state_utility <- q_remission/12
          }
        }
        
        # Apply discounting and add to cycle totals
        cycle_costs[i, cycle_num] <- cycle_costs[i, cycle_num] + 
          apply_discount(state_pop * state_cost, disc_rate, current_year)
        cycle_qalys[i, cycle_num] <- cycle_qalys[i, cycle_num] + 
          apply_discount(state_pop * state_utility, disc_rate, current_year)
        
        # Steroid-related adverse events
        steroid_ae_probs <- get_steroid_ae_probs_by_state(
          cycle = cycle_num,
          state = "Initial remission",
          arm = i,
          duration = steroid_effect_duration_months,
          time_in_state = cycle_num
        )
        
        for (event in all_events) {
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
        
        # State costs and utilities
        state_utility <- (q_relapse - q_decrement_steroid)/12
        
        if (i == 1) {  # Treatment arm
          state_cost <- c_steroid_treat_relapse + c_steroid_monitoring
        } else {  # Control arm
          state_cost <- c_steroid_control_relapse + c_steroid_monitoring
        }
        
        # Apply discounting and add to cycle totals
        cycle_costs[i, cycle_num] <- cycle_costs[i, cycle_num] +
          apply_discount(state_pop * state_cost, disc_rate, current_year)
        cycle_qalys[i, cycle_num] <- cycle_qalys[i, cycle_num] +
          apply_discount(state_pop * state_utility, disc_rate, current_year)
        
        # GCA-related complications (Vision Loss): only occur in relapse state
        n_vl_incident <- pop["Relapse", j-1, i] * p_vl_monthly
        # Update the incident event log by adding to the existing count for this cycle
        event_log$vl[i, cycle_num] <- event_log$vl[i, cycle_num] + n_vl_incident
        
        # Steroid-related adverse events in relapse
        steroid_ae_probs <- get_steroid_ae_probs_by_state(
          cycle = cycle_num,
          state = "Relapse",
          arm = i,
          duration = steroid_effect_duration_months,
          time_in_state = 1
        )
        
        for (event in all_events) {
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
          
          # State costs and utilities
          if (i == 1) {  # Treatment arm
            
            if (tunnel_idx <= 7) {  # steroid tapering applies for first 7 tunnel states (8-month schedule minus 1 relapse cycle) 
              state_utility <- (q_remission - q_decrement_steroid)/12
              state_cost <- c_steroid_monitoring + 
                ifelse(tunnel_idx == 1, c_steroid_treat_sub_remission, 0)
            } else {
              state_utility <- q_remission/12
              state_cost <- 0
            }
            
          } else {  # Control arm
            
            if (tunnel_idx <= 10) {  # steroid tapering applies for first 10 tunnel states (11-month schedule minus 1 relapse cycle)
              state_utility <- (q_remission - q_decrement_steroid)/12
              state_cost <- c_steroid_monitoring + 
                ifelse(tunnel_idx == 1, c_steroid_control_sub_remission, 0)
            } else {
              state_utility <- q_remission/12
              state_cost <- 0
            }
          }
          
          # Apply discounting and add to cycle totals
          cycle_costs[i, cycle_num] <- cycle_costs[i, cycle_num] + 
            apply_discount(state_pop * state_cost, disc_rate, current_year)
          cycle_qalys[i, cycle_num] <- cycle_qalys[i, cycle_num] + 
            apply_discount(state_pop * state_utility, disc_rate, current_year)
          
          # Steroid-related adverse events
          steroid_ae_probs <- get_steroid_ae_probs_by_state(
            cycle = cycle_num,
            state = "Subsequent remission",
            arm = i,
            duration = steroid_effect_duration_months,
            time_in_state = tunnel_idx
          )
          
          for (event in all_events) {
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
        
        # NEW CASES - YEAR 1: Apply year 1 difference for first 12 months
        if (event_log[[event]][i, cycle_num] > 1e-10) {
          new_cases <- event_log[[event]][i, cycle_num]
          
          event_params <- get_event_params(event)
          
          # Apply the DIFFERENCE between year 1 and year 2+ for 12 months
          cost_y1_diff <- event_params$cost_y1 - event_params$cost_y2on
          decrement_y1_diff <- event_params$decrement_y1 - event_params$decrement_y2on
          
          # Apply difference for 12 months starting this cycle
          for (month_in_year1 in 0:11) {
            if (cycle_num + month_in_year1 <= time_horizon_months) {
              future_cycle <- cycle_num + month_in_year1
              
              discounted_cost <- apply_discount(new_cases * cost_y1_diff, 
                                                disc_rate, 
                                                floor((future_cycle - 1)/12) + 1)
              discounted_qaly <- apply_discount(new_cases * decrement_y1_diff, 
                                                disc_rate, 
                                                floor((future_cycle - 1)/12) + 1)
              
              event_costs[i, future_cycle] <- event_costs[i, future_cycle] + discounted_cost 
              
              
              event_qalys[i, future_cycle] <- event_qalys[i, future_cycle] - discounted_qaly 
              
              # Per-event tracking - same values
              idx <- match(event, all_events)
              event_costs_by_event[idx, i, future_cycle] <- 
                event_costs_by_event[idx, i, future_cycle] + discounted_cost
              event_qalys_by_event[idx, i, future_cycle] <- 
                event_qalys_by_event[idx, i, future_cycle] - discounted_qaly
              
            }
          }
        }
        
        # ALL EXISTING CASES - YEAR 2+: Apply year 2+ baseline costs every cycle
        existing_cumulative <- cumulative_populations[[event]][i, j-1]
        
        if (existing_cumulative > 1e-10) {
          event_params <- get_event_params(event)
          
          # Apply year 2+ costs/utilities to all existing cases (monthly basis)
          discounted_cost <- apply_discount(existing_cumulative * event_params$cost_y2on, 
                                            disc_rate, current_year)
          discounted_qaly <- apply_discount(existing_cumulative * event_params$decrement_y2on, 
                                            disc_rate, current_year)
          
          event_costs[i, cycle_num] <- event_costs[i, cycle_num] + discounted_cost
          event_qalys[i, cycle_num] <- event_qalys[i, cycle_num] - discounted_qaly
          
          # Per-event tracking - same values
          idx <- match(event, all_events)
          event_costs_by_event[idx, i, cycle_num] <- 
            event_costs_by_event[idx, i, cycle_num] + discounted_cost
          event_qalys_by_event[idx, i, cycle_num] <- 
            event_qalys_by_event[idx, i, cycle_num] - discounted_qaly
          
        }
      }
      
    }
  }
  
  # Calculate total costs and QALYs
  total_costs <- cycle_costs + event_costs
  total_qalys <- cycle_qalys + event_qalys
  
  # Create summary populations for the original 4 states for reporting (need to update)
  pop_summary <- array(0, 
                       dim = c(4, time_horizon_months + 1, n_arms),
                       dimnames = list(state = c("Initial remission", "Relapse", "Subsequent remission", "Dead"), 
                                       cycle = 0:time_horizon_months, 
                                       arm = arms))
  
  # Aggregate tunnel states back to main states
  for (i in 1:n_arms) {
    for (j in 1:(time_horizon_months + 1)) {
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
  icer <- if (abs(inc_qalys) < 1e-10) NA else (inc_cost / inc_qalys)
  
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

summarize_adverse_events <- function(results, time_horizon_months = 5*12) {
  
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
  alive_pop <- array(0, dim = c(n_arms, time_horizon_months))
  for (i in 1:n_arms) {
    for (j in 1:time_horizon_months) {
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
    GC_alone_cost    = round(total_event_costs_by_event[, "2"], 2),
    Difference_cost  = round(total_event_costs_by_event[, "1"] - total_event_costs_by_event[, "2"], 2),
    TCZ_plus_GC_QALY = round(total_event_qalys_by_event[, "1"], 4),
    GC_alone_QALY    = round(total_event_qalys_by_event[, "2"], 4),
    Difference_QALY  = round(total_event_qalys_by_event[, "1"] - total_event_qalys_by_event[, "2"], 4),
    check.names = FALSE
  )
  
  return(list(
    total_incident_events = total_incident_events_summary,
    event_rates_per_100py = event_rates_per_100py_summary,
    total_patient_years = total_patient_years,
    per_event_cost_qaly = per_event_summary
  ))
  
}
