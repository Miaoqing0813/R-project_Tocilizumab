################
# RUN ANALYSIS
###############
rm(list=ls())

source("01_data/01_model_inputs.R")
source("01_data/02_survival analyses_all.R")

# Depending on whether TCZ continue until planned treatment duration
source("02_decision model/01_Markov_model.R")

# Deterministic CEA results 
############################

subgroup = "itt"  # Options: "newlydiag", "relapsed", "itt"
treatment_duration_months <- 12
treatment_duration_weeks <- 52
time_horizon_months <- 30*12
time_horizon_weeks <- 30*52
steroid_effect_duration_months <- 30*12 # base case: 30 years of steroid effect; use 5 years, 10 years and 20 years in the scenario analyses
dist_beyond <- "weibull" # base case: weibull; use exp for sensitivity analyses
q_decrement_steroid_duration_months <- 2*12 # base case: 2 years of decrement for steroid-related adverse events; use 5 years and 10 years in the scenario analyses
# source("02_decision model/01a_Markov_model_TCZ_relapse.R") # use when continuing TCZ at relapse

# Transition probabilities to first flare 
surv_results <- generate_first_flare_transition_probs(
  subgroup = subgroup,
  treatment_duration_weeks = treatment_duration_weeks,
  time_horizon_weeks = time_horizon_weeks,
  time_horizon_months = time_horizon_months,
  dist_beyond = dist_beyond,
  generate_plots = TRUE,
  verbose = TRUE
)

print(surv_results$trans_prob_treat_plot)
print(surv_results$trans_prob_control_plot)
summary(surv_results$models$control_beyond, type = "survival", t = 260, ci = TRUE)[[1]] # 5-year survival for the comparison group
trans_probs_first_flare_treat <- surv_results$trans_probs_treat
trans_probs_first_flare_control <- surv_results$trans_probs_control

baseline_results <- run_markov_model(use_baseline = TRUE, 
                                     subgroup = subgroup,
                                     time_horizon_weeks = time_horizon_weeks,
                                     time_horizon_months = time_horizon_months,
                                     treatment_duration_months = treatment_duration_months,
                                     steroid_effect_duration_months = steroid_effect_duration_months, 
                                     q_decrement_steroid_duration_months = q_decrement_steroid_duration_months
                                     )

# Print cost effectiveness results
cat("\n=== BASELINE RESULTS (ITT: 1-year treatment) ===\n")
cat(paste("Total cost (Tocilizumab + GC): £", round(baseline_results$total_costs_sum[1], 2), "\n"))
cat(paste("Total QALYs (Tocilizumab + GC):", round(baseline_results$total_qalys_sum[1], 4), "\n"))
cat(paste("Total cost (GC alone): £", round(baseline_results$total_costs_sum[2], 2), "\n"))
cat(paste("Total QALYs (GC alone):", round(baseline_results$total_qalys_sum[2], 4), "\n"))
cat(paste("Incremental cost: £", round(baseline_results$inc_cost, 2), "\n"))
cat(paste("Incremental QALYs:", round(baseline_results$inc_qalys, 4), "\n"))
cat(paste("ICER: £", round(baseline_results$icer, 2), " per QALY\n"))

# Print adverse event summary
baseline_ae_summary <- summarize_adverse_events(baseline_results)

cat("\n=== TOTAL INCIDENT EVENTS ===\n")
cat("Clinical interpretation guidance:\n")
cat("- <1.0: Most patients do not experience this event\n")
cat("- ~1.0: On average, each patient experiences this event once\n")
cat("- >1.0: Some patients experience multiple episodes of this event\n\n")
print(baseline_ae_summary$total_incident_events, row.names = FALSE)

cat("\n=== EVENT RATES PER 100 PATIENT-YEARS ===\n")
print(baseline_ae_summary$event_rates_per_100py, row.names = FALSE)
cat(paste("\nTotal patient-years - TCZ+GC:", round(baseline_ae_summary$total_patient_years[1], 1), 
          "| GC alone:", round(baseline_ae_summary$total_patient_years[2], 1), "\n"))

cat("\n=== LIFETIME PER-EVENT COST AND QALY IMPACTS ===\n")
print(baseline_ae_summary$per_event_cost_qaly, row.names = FALSE)

# Calculate total event costs
total_event_costs <- rowSums(baseline_results$event_costs)
cat("\n=== TOTAL EVENT COSTS ===\n")
cat(paste("TCZ + GC arm: £", round(total_event_costs[1], 0), "\n"))
cat(paste("GC alone arm: £", round(total_event_costs[2], 0), "\n"))
cat(paste("Difference: £", round(total_event_costs[1] - total_event_costs[2], 0), "\n"))
cat(paste("Cost savings from reduced adverse events: £", round(total_event_costs[2] - total_event_costs[1], 0), "\n"))

# Generate PSA parameters
n_sim <- 1000
source("01_data/03_psa_inputs.R")

psa_params_all <- generate_psa_params(n_sim = n_sim, 
                                      seed = 12345, 
                                      time_horizon_months = time_horizon_months)

# Create arrays to store results
psa_total_costs <- matrix(NA, nrow = n_sim, ncol = n_arms)
psa_total_qalys <- matrix(NA, nrow = n_sim, ncol = n_arms)
psa_icer <- numeric(n_sim)

# Store adverse event outcomes from PSA
psa_incident_events <- array(NA, 
                             dim = c(n_sim, 7, n_arms),
                             dimnames = list(sim = NULL, 
                                             event = c("Vision Loss", "Stroke", "Myocardial Infarction", 
                                                       "Heart Failure", "Fractures", "Diabetes", "Infection"),
                                             arm = c("TCZ_plus_GC", "GC_alone")))
psa_event_rates_per_100py <- psa_per_event_costs <- psa_per_event_qalys <- psa_incident_events

# Run PSA simulations
n_cores <- detectCores() - 1  # Leave one core free
cl <- makeCluster(n_cores)

# Source all necessary files on each cluster node
# This automatically makes all functions and variables available
clusterCall(cl, function() {
  source("01_data/01_model_inputs.R")
  source("02_decision model/01_Markov_model.R")
})

# Now you only need to export the iteration-specific data
clusterExport(cl, c("psa_params_all", 
                    "subgroup", 
                    "time_horizon_weeks", 
                    "time_horizon_months", 
                    "treatment_duration_months", 
                    "treatment_duration_weeks",
                    "steroid_effect_duration_months",
                    "q_decrement_steroid_duration_months",
                    "trans_probs_first_flare_treat",
                    "trans_probs_first_flare_control"),
              envir = environment())

# Run PSA in parallel
psa_results_list <- parLapply(cl, 1:n_sim, function(sim) {
  sim_params <- list(
    scalar_params = as.list(psa_params_all$scalar_params[sim, ]),
    vector_params = lapply(psa_params_all$vector_params, function(x) x[sim, ])
  )
  
  results <- run_markov_model(sim_params, 
                              use_baseline = FALSE, 
                              subgroup = subgroup,
                              time_horizon_weeks = time_horizon_weeks,
                              time_horizon_months = time_horizon_months,
                              treatment_duration_months = treatment_duration_months,
                              steroid_effect_duration_months = steroid_effect_duration_months)
  
  ae_summary <- summarize_adverse_events(results, time_horizon_months = time_horizon_months)
  
  return(list(
    total_costs = results$total_costs_sum,
    total_qalys = results$total_qalys_sum,
    icer = results$icer,
    total_incident_events = ae_summary$total_incident_events,
    event_rates_per_100py = ae_summary$event_rates_per_100py,
    per_event_cost_qaly = ae_summary$per_event_cost_qaly,
    total_patient_years = ae_summary$total_patient_years
  ))
})

stopCluster(cl)

# Convert results back to arrays
psa_total_costs <- do.call(rbind, lapply(psa_results_list, function(x) x$total_costs))
psa_total_qalys <- do.call(rbind, lapply(psa_results_list, function(x) x$total_qalys))
psa_icer <- sapply(psa_results_list, function(x) x$icer)

inc_costs_psa <- psa_total_costs[, 1] - psa_total_costs[, 2]
inc_qalys_psa <- psa_total_qalys[, 1] - psa_total_qalys[, 2]
psa_icer_mean <- mean(inc_costs_psa)/mean(inc_qalys_psa)

# Summary statistics
cat("\n=== PSA RESULTS SUMMARY ===\n")
cat(paste("Mean Incremental Cost: £", round(mean(inc_costs_psa), 2), 
          " (95% CI: £", round(quantile(inc_costs_psa, 0.025), 2), " to £", 
          round(quantile(inc_costs_psa, 0.975), 2), ")\n", sep = ""))

cat(paste("Mean Incremental QALYs:", round(mean(inc_qalys_psa), 4), 
          " (95% CI:", round(quantile(inc_qalys_psa, 0.025), 4), "to", 
          round(quantile(inc_qalys_psa, 0.975), 4), ")\n"))

cat(paste("Mean ICER: £", round(psa_icer_mean, 2), 
          " per QALY\n", sep = ""))

# Store adverse event counts
# Populate psa_incident_events from psa_results_list
event_names <- c("Vision Loss", "Stroke", "Myocardial Infarction", 
                 "Heart Failure", "Fractures", "Diabetes", "Infection")

for (sim in 1:n_sim) {
  incident_df <- psa_results_list[[sim]]$total_incident_events
  rates_df <- psa_results_list[[sim]]$event_rates_per_100py
  cost_qaly_df <- psa_results_list[[sim]]$per_event_cost_qaly
  
  for (event_idx in 1:7) {
    event_name <- event_names[event_idx]
    
    # Find matching row in dataframes
    row_idx <- which(incident_df$Event == event_name)
    
    if (length(row_idx) > 0) {
      # Populate psa_incident_events
      psa_incident_events[sim, event_idx, 1] <- incident_df[row_idx, "TCZ_plus_GC"]
      psa_incident_events[sim, event_idx, 2] <- incident_df[row_idx, "GC_alone"]
      
      # Populate psa_event_rates_per_100py
      psa_event_rates_per_100py[sim, event_idx, 1] <- rates_df[row_idx, "TCZ_plus_GC"]
      psa_event_rates_per_100py[sim, event_idx, 2] <- rates_df[row_idx, "GC_alone"]
      
      # Populate psa_per_event_costs
      psa_per_event_costs[sim, event_idx, 1] <- cost_qaly_df[row_idx, "TCZ_plus_GC_cost"]
      psa_per_event_costs[sim, event_idx, 2] <- cost_qaly_df[row_idx, "GC_alone_cost"]
      
      # Populate psa_per_event_qalys
      psa_per_event_qalys[sim, event_idx, 1] <- cost_qaly_df[row_idx, "TCZ_plus_GC_QALY"]
      psa_per_event_qalys[sim, event_idx, 2] <- cost_qaly_df[row_idx, "GC_alone_QALY"]
    }
  }
}

# Now you can use all arrays for summary statistics:
cat("\n=== PSA ADVERSE EVENT SUMMARY ===\n")
for (event_idx in 1:7) {
  event_name <- dimnames(psa_incident_events)$event[event_idx]
  mean_tcz <- mean(psa_incident_events[, event_idx, 1])
  mean_gc <- mean(psa_incident_events[, event_idx, 2])
  ci_tcz <- quantile(psa_incident_events[, event_idx, 1], c(0.025, 0.975))
  ci_gc <- quantile(psa_incident_events[, event_idx, 2], c(0.025, 0.975))
  cat(paste0(event_name, ":\n"))
  cat(paste0("  TCZ+GC: ", round(mean_tcz, 3), " (95% CI: ",
             round(ci_tcz[1], 3), " to ", round(ci_tcz[2], 3), ")\n"))
  cat(paste0("  GC alone: ", round(mean_gc, 3), " (95% CI: ",
             round(ci_gc[1], 3), " to ", round(ci_gc[2], 3), ")\n"))
  cat(paste0("  Difference: ", round(mean_tcz - mean_gc, 3), "\n\n"))
}

cat("\n=== PSA EVENT RATES (per 100 patient-years) ===\n")
for (event_idx in 1:7) {
  event_name <- dimnames(psa_event_rates_per_100py)$event[event_idx]
  mean_tcz <- mean(psa_event_rates_per_100py[, event_idx, 1])
  mean_gc <- mean(psa_event_rates_per_100py[, event_idx, 2])
  ci_tcz <- quantile(psa_event_rates_per_100py[, event_idx, 1], c(0.025, 0.975), na.rm = TRUE)
  ci_gc <- quantile(psa_event_rates_per_100py[, event_idx, 2], c(0.025, 0.975), na.rm = TRUE)
  cat(paste0(event_name, ":\n"))
  cat(paste0("  TCZ+GC: ", round(mean_tcz, 2), " (95% CI: ",
             round(ci_tcz[1], 2), " to ", round(ci_tcz[2], 2), ")\n"))
  cat(paste0("  GC alone: ", round(mean_gc, 2), " (95% CI: ",
             round(ci_gc[1], 2), " to ", round(ci_gc[2], 2), ")\n"))
  cat(paste0("  Difference: ", round(mean_tcz - mean_gc, 2), "\n\n"))
}

cat("\n=== PSA PER-EVENT COSTS ===\n")
for (event_idx in 1:7) {
  event_name <- dimnames(psa_per_event_costs)$event[event_idx]
  mean_tcz <- mean(psa_per_event_costs[, event_idx, 1])
  mean_gc <- mean(psa_per_event_costs[, event_idx, 2])
  ci_tcz <- quantile(psa_per_event_costs[, event_idx, 1], c(0.025, 0.975), na.rm = TRUE)
  ci_gc <- quantile(psa_per_event_costs[, event_idx, 2], c(0.025, 0.975), na.rm = TRUE)
  cat(paste0(event_name, ":\n"))
  cat(paste0("  TCZ+GC: £", round(mean_tcz, 2), " (95% CI: £",
             round(ci_tcz[1], 2), " to £", round(ci_tcz[2], 2), ")\n"))
  cat(paste0("  GC alone: £", round(mean_gc, 2), " (95% CI: £",
             round(ci_gc[1], 2), " to £", round(ci_gc[2], 2), ")\n"))
  cat(paste0("  Difference: £", round(mean_tcz - mean_gc, 2), "\n\n"))
}

cat("\n=== PSA PER-EVENT QALYs ===\n")
for (event_idx in 1:7) {
  event_name <- dimnames(psa_per_event_qalys)$event[event_idx]
  mean_tcz <- mean(psa_per_event_qalys[, event_idx, 1])
  mean_gc <- mean(psa_per_event_qalys[, event_idx, 2])
  ci_tcz <- quantile(psa_per_event_qalys[, event_idx, 1], c(0.025, 0.975), na.rm = TRUE)
  ci_gc <- quantile(psa_per_event_qalys[, event_idx, 2], c(0.025, 0.975), na.rm = TRUE)
  cat(paste0(event_name, ":\n"))
  cat(paste0("  TCZ+GC: ", round(mean_tcz, 4), " (95% CI: ",
             round(ci_tcz[1], 4), " to ", round(ci_tcz[2], 4), ")\n"))
  cat(paste0("  GC alone: ", round(mean_gc, 4), " (95% CI: ",
             round(ci_gc[1], 4), " to ", round(ci_gc[2], 4), ")\n"))
  cat(paste0("  Difference: ", round(mean_tcz - mean_gc, 4), "\n\n"))
}

# Cost-effectiveness plane
png("03_results/CE_Plane_Plot.png", width = 2000, height = 1600, res = 300)

max_x <- max(abs(inc_qalys_psa), na.rm = TRUE)
max_y <- max(abs(inc_costs_psa), na.rm = TRUE)

plot(inc_qalys_psa, inc_costs_psa, 
     xlab = "Incremental QALYs", 
     ylab = "Incremental Costs (£)",
     main = "Cost-Effectiveness Plane",
     pch = 20, col = rgb(0, 0, 1, 0.3), 
     xlim = c(-max_x, max_x), 
     ylim = c(-max_y, max_y)) 
abline(h = 0, v = 0, lty = 2)

wtp_thresholds <- c(25000, 35000)
for (wtp in wtp_thresholds) {
  abline(a = 0, b = wtp, col = "red", lty = 2)
}
legend("topleft", legend = c("£25,000/QALY", "£35,000/QALY"), 
       col = "red", lty = 2)

dev.off()

# CEAC
png("03_results/CEAC_Plot.png", width = 2000, height = 1600, res = 300)

wtp_range <- seq(0, 50000, by = 1000)
prob_ce <- sapply(wtp_range, function(wtp) {
  mean((inc_costs_psa / inc_qalys_psa) <= wtp)
})

plot(wtp_range, prob_ce, type = "l", lwd = 2,
     xlab = "Willingness to Pay (£/QALY)",
     ylab = "Probability Cost-Effective",
     main = "Cost-Effectiveness Acceptability Curve")
abline(h = 0.5, lty = 2, col = "gray")
abline(v = c(25000, 35000), lty = 2, col = "red")

dev.off()

# Save results in Excel
source("02_decision model/03_export_Excel.R")


