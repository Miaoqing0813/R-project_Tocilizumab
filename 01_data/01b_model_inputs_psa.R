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

generate_psa_params <- function(n_sim = 10, seed=12345) {
  
  set.seed(seed)
  
  # List of all cost variable names (log-normal distribution)
  cost_vars <- c("c_TCZ", "c_TCZ_PAS", "c_steroid_control_init_remission", 
                 "c_steroid_treat_init_remission", "c_steroid_control_relapse", 
                 "c_steroid_treat_relapse", "c_steroid_control_sub_remission", 
                 "c_steroid_treat_sub_remission", "c_steroid_monitoring", 
                 "c_TCZ_monitoring",
                 "c_vl_y1", "c_vl_y2on", "c_stroke_y1", "c_mi_y1", "c_hf_y1",
                 "c_stroke_y2on", "c_mi_y2on", "c_hf_y2on", "c_fracture_y1",
                 "c_fracture_y2on", "c_bp", "c_dm", "c_infection")
  
  # List of all utility variable names (beta distribution)
  utility_vars <- c("q_remission", "q_relapse", "q_vl", "q_stroke_y1", 
                    "q_mi_y1", "q_hf_y1", "q_stroke_y2on", "q_mi_y2on", 
                    "q_hf_y2on", "q_fracture_y1", "q_fracture_y2on", "q_dm")
  
  # List of all probability variable names (beta distribution)
  prob_vars <- c("p_vl", "p_stroke", "p_mi", "p_hf", "p_fractures", 
                 "p_dm", "p_infection")
  
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
  
  # Return both scalar and vector parameters as a list
  return(list(
    scalar_params = l_params_psa,
    vector_params = l_params_psa_vectors
  ))
}



