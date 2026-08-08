#####################
# Model validation # 
#####################

# Quick checks 
################

# All input parameters are used in the analysis (i.e. no redundant parameters are retained in the model)

# Check transition probabilities: after running get_tpm(), verify each row sums to 1
row_sums <- rowSums(tpm)
all(abs(row_sums - 1) < 1e-10)  # Should be TRUE

# Check whether total population distribution sum to one
for (i in 1:n_arms) {
  total_pop_per_cycle <- colSums(baseline_results$pop[, , i])
  print(paste("Arm", i, "starts with", total_pop_per_cycle[1], 
              "and ends with", total_pop_per_cycle[length(total_pop_per_cycle)]))
}

# Check there's no negative cost
baseline_results[["cycle_costs"]]
baseline_results[["event_costs"]]

# Compare total QALYs and total life years
baseline_results$total_qalys_sum
sum(baseline_results$pop_summary[c("Initial remission", "Relapse", "Subsequent remission"), 
                                 2:ncol(baseline_results$pop_summary), 1]) / 12
sum(baseline_results$pop_summary[c("Initial remission", "Relapse", "Subsequent remission"), 
                                 2:ncol(baseline_results$pop_summary), 2]) / 12

# Proportion of patients in the death state should not decrease over time
baseline_results$pop_summary["Dead", , 1]
death_diff_arm1 <- diff(baseline_results$pop_summary["Dead", , 1])
any(death_diff_arm1 < 0)

baseline_results$pop_summary["Dead", , 2]
death_diff_arm1 <- diff(baseline_results$pop_summary["Dead", , 2])
any(death_diff_arm1 < 0)

# ============================
# PREDICTIONS OF AES VALIDATION
# =============================

## The Per-100 Patient-Year Rate Test
## Cumulative sum of new diabetes cases up to month 12 for both arms
model_1yr_risk_arm1 <- sum(baseline_results$event_log$dm[1, 1:12]) * 100
model_1yr_risk_arm2 <- sum(baseline_results$event_log$dm[2, 1:12]) * 100
cat("Model Predicted 1-Year Diabetes Risk (Arm 1):", round(model_1yr_risk_arm1, 2), "%\n")
cat("Model Predicted 1-Year Diabetes Risk (Arm 2):", round(model_1yr_risk_arm2, 2), "%\n")

# The sum of all CEAC values at each WTP threshold should equal 1
wtp_thresholds <- seq(0, 50000, by = 1000)
validation_matrix <- data.frame(WTP = wtp_thresholds, Arm1_Prob = NA, Arm2_Prob = NA, Total = NA)
  
for (i in seq_along(wtp_thresholds)) {
  lambda <- wtp_thresholds[i]
  
  # Calculate Net Monetary Benefit for each iteration in the PSA
  nmb_arm1 <- (psa_total_qalys[, 1] * lambda) - psa_total_costs[, 1]
  nmb_arm2 <- (psa_total_qalys[, 2] * lambda) - psa_total_costs[, 2]
  
  # Count how often each arm wins at this specific threshold
  validation_matrix$Arm1_Prob[i] <- mean(nmb_arm1 > nmb_arm2)
  validation_matrix$Arm2_Prob[i] <- mean(nmb_arm2 >= nmb_arm1)
  
  # HORIZONTAL SUMMATION CHECK
  validation_matrix$Total[i] <- validation_matrix$Arm1_Prob[i] + validation_matrix$Arm2_Prob[i]
}

# Print the table to verify
print(validation_matrix) # the total should be one!
  
# ======================
# EVENT COST VALIDATION
# ======================
cat("\n=== CORRECTED EVENT COST VALIDATION ===\n\n")

# For each event: Check cost per incident (proxy for reasonableness)
for (idx in 1:length(events)) {
  event <- events[idx]
  label <- event_labels[idx]
  
  n_tcz <- incident_events[idx, "TCZ_plus_GC"]
  n_gc <- incident_events[idx, "GC_alone"]
  cost_tcz <- per_event_costs[idx, "TCZ_plus_GC_cost"]
  cost_gc <- per_event_costs[idx, "GC_alone_cost"]
  
  # Cost per incident (proxy metric)
  cost_per_incident_tcz <- if (n_tcz > 0) cost_tcz / n_tcz else 0
  cost_per_incident_gc <- if (n_gc > 0) cost_gc / n_gc else 0
  
  cat(paste(label, ":\n"))
  cat(paste("  TCZ: ", round(n_tcz, 4), " incidents @ £", 
            round(cost_per_incident_tcz, 0), "/incident = £", 
            round(cost_tcz, 0), "\n"))
  cat(paste("  GC:  ", round(n_gc, 4), " incidents @ £", 
            round(cost_per_incident_gc, 0), "/incident = £", 
            round(cost_gc, 0), "\n"))
  
  # Sanity check: GC cost per incident should be similar to TCZ
  # (same event, just different frequencies due to steroid exposure)
  ratio <- cost_per_incident_gc / cost_per_incident_tcz
  
  if (ratio < 0.7 | ratio > 1.4) {
    cat(paste("  ⚠️ Cost per incident ratio (GC/TCZ): ", round(ratio, 2), 
              " - CHECK CALCULATION\n"))
  } else {
    cat(paste("  ✓ Cost per incident ratio (GC/TCZ): ", round(ratio, 2), "\n"))
  }
  cat("\n")
}
  
# Compare model results with back of the envelope calculations
  
# The PSA results should be fairly similar to the deterministic results
  ## Compare baseline with PSA sampled means
  cost_vars <- c("c_steroid_monitoring","c_TCZ_monitoring",
                 "c_vl_y1","c_vl_y2on","c_stroke_y1","c_mi_y1","c_hf_y1",
                 "c_stroke_y2on","c_mi_y2on","c_hf_y2on","c_fracture_y1",
                 "c_fracture_y2on","c_dm","c_infection")

  for (v in cost_vars) {
    baseline_val <- get(v)
    sampled_mean <- mean(psa_params_all$scalar_params[[v]])
    cat(v, " baseline:", baseline_val, " sampled_mean:", sampled_mean, "\n")
  }
  
  utility_vars <- c("q_remission", "q_relapse", "q_decrement_dm", 
                    "q_multiplier_vl", "q_multiplier_stroke_y1", "q_multiplier_mi_y1", 
                    "q_multiplier_hf_y1", "q_multiplier_stroke_y2on", "q_multiplier_mi_y2on", 
                    "q_multiplier_hf_y2on", "q_multiplier_fracture_y1", "q_multiplier_fracture_y2on")
  
  for (v in utility_vars) {
    baseline_val <- get(v)
    sampled_mean <- mean(psa_params_all$scalar_params[[v]])
    cat(v, " baseline:", baseline_val, " sampled_mean:", sampled_mean, "\n")
  }
  
  # Check transition probs (first-flare)
  print(data.frame(month = 1:12, 
                   baseline_tp = trans_probs_first_flare_treat$trans_prob[1:12], 
                   sampled_tp_mean = colMeans(psa_params_all$vector_params$trans_probs_first_flare_treat)[1:12]))
  cat(" baseline:", probs_subseq_flare_itt_treat, 
      " sampled_mean:", mean(psa_params_all$scalar_params$probs_subseq_flare_itt_treat))
  cat(" baseline:", probs_subseq_flare_newlydiag_treat, 
      " sampled_mean:", mean(psa_params_all$scalar_params$probs_subseq_flare_newlydiag_treat))
  cat(" baseline:", probs_subseq_flare_relapsed_treat, 
      " sampled_mean:", mean(psa_params_all$scalar_params$probs_subseq_flare_relapsed_treat))
  
  ## Reproduce baseline by passing PSA means back into the model
  ## Should get almost the same totals when using PSA means
  scalar_means <- colMeans(psa_params_all$scalar_params)
  
  vector_means <- lapply(psa_params_all$vector_params, function(mat) {
    if (is.matrix(mat) || is.data.frame(mat)) {
      as.numeric(colMeans(mat))
    } else {
      as.numeric(mat)  # if already vector
    }
  })
  
  params_means <- list(scalar_params = as.list(scalar_means),
                       vector_params = vector_means)
  
  res_psa_means <- run_markov_model(params = params_means, use_baseline = FALSE,
                                    subgroup = subgroup,
                                    time_horizon_months = time_horizon_months,
                                    time_horizon_weeks = time_horizon_weeks,
                                    treatment_duration_months = treatment_duration_months,
                                    steroid_effect_duration_months = steroid_effect_duration_months)
  
  print(baseline_results$total_costs)
  print(res_psa_means$total_costs_sum)
  print(baseline_results$total_qalys)
  print(res_psa_means$total_qalys_sum)
  
  # Check the distributin of ICERs
  mean(inc_costs_psa) / mean(inc_qalys_psa) # Ratio of means
  mean(psa_icer[is.finite(psa_icer)], na.rm = TRUE) # Mean of ICERs 
  median(psa_icer[is.finite(psa_icer)], na.rm = TRUE) # Median of ICERs 
  
  summary(inc_qalys_psa)
  table(sign(inc_qalys_psa))           # how many positive / negative / zero-ish
  sum(abs(inc_qalys_psa) < 1e-6)       # count near-zero denominators
  summary(psa_icer)
  quantile(psa_icer[is.finite(psa_icer)], probs = c(0.01,0.1,0.25,0.5,0.75,0.9,0.99))
  # identify extreme sims
  which(abs(psa_icer) > 1e6)           # sims producing huge ICERs
  
  lambda <- 25000
  nmb <- lambda * inc_qalys_psa - inc_costs_psa
  mean_nmb <- mean(nmb)
  prob_cost_effective <- mean(nmb > 0)   # probability TCZ is cost-effective at this WTP
  cat("mean NMB:", mean_nmb, "prob CE:", prob_cost_effective, "\n")
  
  

# Run univariate sensitivity analyses
#####################################
# Increase and decrease the time horizon

# Increase and decrease discount rate for costs and outcomes 
disc_rate <- 0

# Set all cost inputs to zero (or increase/decrease cost inputs)
all_cost_vars <- ls(pattern = "^c_", envir = .GlobalEnv)
for (var_name in all_cost_vars) {
  assign(var_name, 0, envir = .GlobalEnv)
}
cat("Successfully zeroed out", length(all_cost_vars), "cost variables.\n")

# Set all utility values to 1 and AE disutilities to 0
disc_rate <- 0 # don't forget to set discount to be zero

all_utility_decrement <- ls(pattern = "^q_decrement", envir = .GlobalEnv)
for (var_name in all_utility_decrement) {
  assign(var_name, 0, envir = .GlobalEnv)
}

all_utility_multiplier <- ls(pattern = "^q_multiplier", envir = .GlobalEnv)
for (var_name in all_utility_multiplier) {
  assign(var_name, 1, envir = .GlobalEnv)
}

  # Run markov model step by step
  v.q_remission_adjusted <- rep(1, time_horizon_months)
  v.q_relapse_adjusted   <- rep(1, time_horizon_months)
  
  # Total QALYs should equal to total life-years
  total_qalys_sum
  sum(pop_summary[c("Initial remission", "Relapse", "Subsequent remission"), 2:ncol(pop_summary), 1]) / 12
  sum(pop_summary[c("Initial remission", "Relapse", "Subsequent remission"), 2:ncol(pop_summary), 2]) / 12
  
# Set treatment efficacy and AEs the same for both arms
trans_probs_first_flare_treat <- trans_probs_first_flare_control
probs_subseq_flare_newlydiag_treat <- probs_subseq_flare_newlydiag_control
probs_subseq_flare_relapsed_treat  <- probs_subseq_flare_relapsed_control
probs_subseq_flare_itt_treat       <- probs_subseq_flare_itt_control

  # Run markov model step by step  
  # Set steroid-related AE probabilities to zero for both arms
  get_steroid_ae_probs_by_state <- function(cycle, state, arm, 
                                            time_in_state = 1,
                                            duration = steroid_effect_duration_months) {
    
    steroid_ae_probs <- list(
      stroke     = 0,
      mi         = 0,
      hf         = 0,
      dm         = 0,
      infection  = 0,
      fractures  = 0,
      mortality  = v.allcause.mort.monthly[cycle] # Pure background mortality, NO HR multiplier
    )
    
    return(steroid_ae_probs)
  }
  
  # Results check: should be identical QALYs between arms
  total_qalys_sum

# Set mortality to zero
v.allcause.mort.monthly <- rep(0, time_horizon_months)

get_steroid_ae_probs_by_state <- function(cycle, state, arm, 
                                          time_in_state = 1,
                                          duration = steroid_effect_duration_months) {
  
  steroid_ae_probs <- list(
    stroke     = 0,
    mi         = 0,
    hf         = 0,
    dm         = 0,
    infection  = 0,
    fractures  = 0,
    mortality  = 0 # FORCE STEROID-RELATED EXCESS MORTALITY TO ZERO
  )
  
  return(steroid_ae_probs)
}

  # Results check: should be zero for both arms
  print(pop_summary["Dead", ncol(pop_summary), 1])
  print(pop_summary["Dead", ncol(pop_summary), 2])
  
# Set steroid-related AE probabilities to zero for both arms
get_steroid_ae_probs_by_state <- function(cycle, state, arm, 
                                          time_in_state = 1,
                                          duration = steroid_effect_duration_months) {
  
  steroid_ae_probs <- list(
    stroke     = 0,
    mi         = 0,
    hf         = 0,
    dm         = 0,
    infection  = 0,
    fractures  = 0,
    mortality  = v.allcause.mort.monthly[cycle] # Pure background mortality, NO HR multiplier
  )
  
  return(steroid_ae_probs)
}

   # Results check: 
   total_costs_sum # total costs should decrease
   total_qalys_sum # total QALYs should increase
   




