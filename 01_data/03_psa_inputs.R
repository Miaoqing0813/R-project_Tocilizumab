###################################################################################
#### Functions to generate input parameters for probabilistic sensitivity analyses
###################################################################################

# Transform values for mean and standard deviation into parameters 
betaPar <- function(m, s) {
  # m:  Mean of the Beta distribution
  # m: Standard deviation of the Beta distribution
  
  var <- s ^ 2
  alpha <- ((1 - m) / var - 1 / m) * m ^ 2
  beta <- alpha * (1 / m - 1)
  
  return(
    list(alpha = alpha, beta = beta)
  )
}

lognPar <- function(m,s) {
  # m: Mean of Log-Normal distribution
  # s: Standard deiviation of Log-Normal distribution
  
  var <- s^2
  meanlog <- log(m) - 0.5 * log(1 + var/m^2)
  varlog <- log(1 + (var/m^2))
  sdlog <- sqrt(varlog)
  
  return(
    list(meanlog = meanlog, sdlog = sdlog)
  )
}

# Calculate SE from 95% CI bounds
calc_se_from_ci <- function(lower, upper) {
  # Assumes normal approximation: 95% CI = mean ± 1.96*SE
  se <- (upper - lower) / (2 * 1.96)
  return(se)
}

# Function to sample transition probabilities from survival analyses
sample_trans_probs <- function(trans_prob_df, n_sim, time_horizon_months = 30*12) {
  
  # Pre-allocate matrix
  sampled_probs <- matrix(NA, nrow = n_sim, ncol = time_horizon_months)
  
  # Draw a single uniform percentile (0 to 1) per simulation run, rather than independently sampling the probabilities
  # This preserves perfect temporal correlation across the entire survival curve
  u_draws <- runif(n_sim, min = 0, max = 1)
  
  for (cycle in 1:time_horizon_months) {
    mean_prob <- trans_prob_df$trans_prob[cycle]
    lower_prob <- trans_prob_df$trans_prob_lower[cycle]
    upper_prob <- trans_prob_df$trans_prob_upper[cycle]
    
    # Calculate SE from confidence intervals
    se_prob <- calc_se_from_ci(lower_prob, upper_prob)
    
    # Handle edge cases where probability is very close to 0 or 1
    if (mean_prob < 1e-6) {
      # If probability is essentially zero, just use zero
      sampled_probs[, cycle] <- rep(0, n_sim)
    } else if (mean_prob > (1 - 1e-6)) {
      # If probability is essentially one, just use one
      sampled_probs[, cycle] <- rep(1, n_sim)
    } else if (se_prob < 1e-10) {
      # If SE is essentially zero (very tight CI), use the mean
      sampled_probs[, cycle] <- rep(mean_prob, n_sim)
    } else {
      # Get beta distribution parameters for this specific cycle
      beta_params <- betaPar(m = mean_prob, s = se_prob)
      
      # Use Quantile Mapping (Inverse Transform Sampling)
      # Maps our uniform draws to this cycle's exact Beta distribution curve.
      sampled_probs[, cycle] <- qbeta(p = u_draws, 
                                      shape1 = beta_params$alpha, 
                                      shape2 = beta_params$beta)
    }
  }
  
  return(sampled_probs)
}


generate_psa_params <- function(n_sim = 10, seed=12345, time_horizon_months = 30*12) {
  
  set.seed(seed)
  
  # List of all cost variable names (log-normal distribution)
  cost_vars <- c("c_steroid_monitoring", "c_TCZ_monitoring",
                 "c_vl_y1", "c_vl_y2on", "c_stroke_y1", "c_mi_y1", "c_hf_y1",
                 "c_stroke_y2on", "c_mi_y2on", "c_hf_y2on", "c_fracture_y1",
                 "c_fracture_y2on", "c_dm", "c_infection")
  
  # List of all utility variable names (beta distribution)
  utility_vars <- c("q_remission", "q_relapse", "q_decrement_dm", "q_decrement_infection",
                    "q_multiplier_vl", "q_multiplier_stroke_y1", "q_multiplier_mi_y1", 
                    "q_multiplier_hf_y1", "q_multiplier_stroke_y2on", "q_multiplier_mi_y2on", 
                    "q_multiplier_hf_y2on", "q_multiplier_fracture_y1", "q_multiplier_fracture_y2on")
  
  # List of all probability variable names (beta distribution)
  prob_vars <- c("p_vl", "p_stroke", "p_mi", "p_hf", "p_fractures", "p_dm", "p_infection", 
                 "probs_subseq_flare_newlydiag_control", "probs_subseq_flare_newlydiag_treat",
                 "probs_subseq_flare_relapsed_control", "probs_subseq_flare_relapsed_treat", 
                 "probs_subseq_flare_itt_control", "probs_subseq_flare_itt_treat")
  
  # Get vector variable names
  all_vars_env <- ls(pattern = "^v\\.steroid\\.hr\\.", envir = parent.frame())
  vector_names <- all_vars_env[!grepl("\\.(lci|uci)$", all_vars_env)]
  
  # Combine all scalar variables
  all_vars <- c(cost_vars, utility_vars, prob_vars)
  
  # Create a data frame to store PSA parameters for scalar variables
  l_params_psa <- data.frame(matrix(nrow = n_sim, ncol = length(all_vars)))
  colnames(l_params_psa) <- all_vars
  
  # Sample values for each cost variable (log-normal)
  for (var_name in cost_vars) {
    mean_val <- get(var_name, envir = parent.frame())
    se_val <- get(paste0(var_name, "_se"), envir = parent.frame())
    
    params <- lognPar(m = mean_val, s = se_val)
    
    l_params_psa[[var_name]] <- rlnorm(n = n_sim, 
                                       meanlog = params$meanlog, 
                                       sdlog = params$sdlog)
  }
  
  # Sample values for each utility variable (beta)
  for (var_name in utility_vars) {
    mean_val <- get(var_name, envir = parent.frame())
    se_val <- get(paste0(var_name, "_se"), envir = parent.frame())
    
    params <- betaPar(m = mean_val, s = se_val)
    
    l_params_psa[[var_name]] <- rbeta(n = n_sim, 
                                      shape1 = params$alpha, 
                                      shape2 = params$beta)
  }
  
  # Sample values for each probability variable (beta)
  for (var_name in prob_vars) {
    mean_val <- get(var_name, envir = parent.frame())
    se_val <- get(paste0(var_name, "_se"), envir = parent.frame())
    
    params <- betaPar(m = mean_val, s = se_val)
    
    l_params_psa[[var_name]] <- rbeta(n = n_sim, 
                                      shape1 = params$alpha, 
                                      shape2 = params$beta)
  }
  
  # Create a list to store vector-based PSA parameters
  l_params_psa_vectors <- list()
  
  # Vector-based parameters for hazard ratios based on steroid dosages
  for (vec in vector_names) {
    values <- get(vec, envir = parent.frame())
    lci <- get(paste0(vec, ".lci"), envir = parent.frame())
    uci <- get(paste0(vec, ".uci"), envir = parent.frame())
    
    # Use logical indexing to handle changes in vector values
    changes <- c(TRUE, diff(values) != 0)  # TRUE for first row or when values change
    
    # Pre-allocate matrix for PSA simulations for vectors
    vec_psa <- matrix(NA, nrow = n_sim, ncol = length(values))
    
    # Simulate values only for unique rows
    unique_indices <- which(changes) 
    vec_psa[, unique_indices] <- sapply(unique_indices, function(i) {
      exp(rnorm(n_sim, mean = log(values[i]), sd = (log(uci[i]) - log(lci[i])) / 3.92))
    })
    
    # Fill non-unique rows with values from the last unique simulation
    for (j in seq_along(values)) {
      if (!changes[j]) {
        vec_psa[, j] <- vec_psa[, j - 1]
      }
    }
    
    # Store the matrix in the list
    l_params_psa_vectors[[vec]] <- vec_psa
  }
  
  # Treatment arm 
  trans_probs_first_flare_treat_df <- get("trans_probs_first_flare_treat", 
                                                  envir = parent.frame())
  l_params_psa_vectors[["trans_probs_first_flare_treat"]] <- 
    sample_trans_probs(trans_probs_first_flare_treat_df, n_sim, time_horizon_months = time_horizon_months)
  
  # Control arm 
  trans_probs_first_flare_control_df <- get("trans_probs_first_flare_control", 
                                                    envir = parent.frame())
  l_params_psa_vectors[["trans_probs_first_flare_control"]] <- 
    sample_trans_probs(trans_probs_first_flare_control_df, n_sim, time_horizon_months = time_horizon_months)
  
  # Return both scalar and vector parameters as a list
  return(list(
    scalar_params = l_params_psa,
    vector_params = l_params_psa_vectors
  ))
}



