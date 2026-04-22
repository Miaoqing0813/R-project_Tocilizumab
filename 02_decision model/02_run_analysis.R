################
# RUN ANALYSIS
###############
rm(list=ls())

source("01_data/01_model_inputs.R")
source("01_data/02_survival analyses_all.R")
source("02_decision model/01_Markov_model.R")

# Deterministic CEA results 
############################

subgroup = "itt"  # Options: "newlydiag", "relapsed", "itt"
treatment_duration_months <- 12
treatment_duration_weeks <- 52
time_horizon_months <- 5*12
time_horizon_weeks <- 5*52
steroid_effect_duration_months <- 5*12
dist_beyond <- "exp"

# Transition probabilities to first flare 
surv_results_itt <- generate_first_flare_transition_probs(
  subgroup = subgroup,
  treatment_duration_weeks = treatment_duration_weeks,
  time_horizon_weeks = time_horizon_weeks,
  time_horizon_months = time_horizon_months,
  dist_beyond = dist_beyond,
  generate_plots = TRUE,
  verbose = TRUE
)

print(surv_results_itt$trans_prob_treat_plot)
print(surv_results_itt$trans_prob_control_plot)
trans_probs_first_flare_treat <- surv_results_itt$trans_probs_treat
trans_probs_first_flare_control <- surv_results_itt$trans_probs_control

baseline_results <- run_markov_model(use_baseline = TRUE, 
                                     subgroup = subgroup,
                                     time_horizon_weeks = time_horizon_weeks,
                                     time_horizon_months = time_horizon_months,
                                     treatment_duration_months = treatment_duration_months,
                                     steroid_effect_duration_months = steroid_effect_duration_months
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
for (sim in 1:n_sim) {
  
  # Extract parameters for this simulation
  sim_params <- list(
    scalar_params = as.list(psa_params_all$scalar_params[sim, ]),
    vector_params = lapply(psa_params_all$vector_params, function(x) x[sim, ])
  )
  
  # Run model with these parameters
  results <- run_markov_model(sim_params, 
                              use_baseline = FALSE, 
                              subgroup = subgroup,
                              time_horizon_weeks = time_horizon_weeks,
                              time_horizon_months = time_horizon_months,
                              treatment_duration_months = treatment_duration_months,
                              steroid_effect_duration_months = steroid_effect_duration_months
  )
  
  # Store results
  psa_total_costs[sim, ] <- results$total_costs_sum
  psa_total_qalys[sim, ] <- results$total_qalys_sum
  psa_icer[sim] <- results$icer
  
  # Store adverse event counts
  ae_summary <- summarize_adverse_events(results)
  psa_incident_events[sim, , 1] <- ae_summary$total_incident_events$TCZ_plus_GC
  psa_incident_events[sim, , 2] <- ae_summary$total_incident_events$GC_alone
  psa_event_rates_per_100py[sim, , 1] <- ae_summary$event_rates_per_100py$TCZ_plus_GC
  psa_event_rates_per_100py[sim, , 2] <- ae_summary$event_rates_per_100py$GC_alone
  psa_per_event_costs[sim, , 1] <- ae_summary$per_event_cost_qaly$TCZ_plus_GC_cost
  psa_per_event_costs[sim, , 2] <- ae_summary$per_event_cost_qaly$GC_alone_cost
  psa_per_event_qalys[sim, , 1] <- ae_summary$per_event_cost_qaly$TCZ_plus_GC_QALY
  psa_per_event_qalys[sim, , 2] <- ae_summary$per_event_cost_qaly$GC_alone_QALY
  
  # Progress indicator
  if (sim %% 100 == 0) {
    cat(paste("Completed", sim, "of", n_sim, "simulations\n"))
  }
}

# Calculate incremental costs and QALYs
inc_costs_psa <- psa_total_costs[, 1] - psa_total_costs[, 2]
inc_qalys_psa <- psa_total_qalys[, 1] - psa_total_qalys[, 2]

# Summary statistics
cat("\n=== PSA RESULTS SUMMARY ===\n")
cat(paste("Mean Incremental Cost: £", round(mean(inc_costs_psa), 2), 
          " (95% CI: £", round(quantile(inc_costs_psa, 0.025), 2), " to £", 
          round(quantile(inc_costs_psa, 0.975), 2), ")\n", sep = ""))

cat(paste("Mean Incremental QALYs:", round(mean(inc_qalys_psa), 4), 
          " (95% CI:", round(quantile(inc_qalys_psa, 0.025), 4), "to", 
          round(quantile(inc_qalys_psa, 0.975), 4), ")\n"))

cat(paste("Mean ICER: £", round(mean(psa_icer[is.finite(psa_icer)]), 2), 
          " per QALY\n", sep = ""))

# Summary of adverse events from PSA
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

wtp_thresholds <- c(20000, 30000)
for (wtp in wtp_thresholds) {
  abline(a = 0, b = wtp, col = "red", lty = 2)
}
legend("topleft", legend = c("£20,000/QALY", "£30,000/QALY"), 
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
abline(v = c(20000, 30000), lty = 2, col = "red")

dev.off()

# Save results in Excel
source("02_decision model/03_export_excel.R")





# Run all scenarios together
##############################

scenarios <- data.frame(
  scenario_name = c("Base case: Newly diagnosed, 52 weeks",
                    "Scenario 1: Newly diagnosed, 104 weeks",
                    "Scenario 2: Relapsed, 52 weeks",
                    "Scenario 3: Relapsed, 104 weeks",
                    "Scenario 4: ITT, 52 weeks", 
                    "Scenario 5: ITT, 104 weeks"),
  subgroup = c("newlydiag", "newlydiag", "relapsed", "relapsed", "itt", "itt"),
  treatment_duration_weeks = c(52, 104, 52, 104, 52, 104),
  treatment_duration_months = c(12, 24, 12, 24, 12, 24),
  stringsAsFactors = FALSE
)

scenario_results <- list()

for (i in 1:nrow(scenarios)) {
  cat("\n========================================\n")
  cat("Running:", scenarios$scenario_name[i], "\n")
  cat("========================================\n")
  
  # Generate transition probabilities for this scenario
  trans_probs <- generate_first_flare_transition_probs(
    subgroup = scenarios$subgroup[i],
    treatment_duration_weeks = scenarios$treatment_duration_weeks[i],
    generate_plots = FALSE,
    verbose = FALSE
  )
  
  # Assign to global environment for Markov model to use
  trans_probs_first_flare_treat <<- trans_probs$trans_probs_treat
  trans_probs_first_flare_control <<- trans_probs$trans_probs_control
  
  # Run Markov model
  markov_results <- run_markov_model(use_baseline = TRUE, 
                                     subgroup = scenarios$subgroup[i],
                                     treatment_duration_months = scenarios$treatment_duration_months[i])
  
  # Generate adverse event summary
  ae_summary <- summarize_adverse_events(markov_results)
  
  # Calculate total event costs
  total_event_costs <- rowSums(markov_results$event_costs)
  
  # Store results
  scenario_results[[scenarios$scenario_name[i]]] <- list(
    markov_results = markov_results,
    trans_probs = trans_probs,
    ae_summary = ae_summary,
    total_event_costs = total_event_costs,
    scenario_info = scenarios[i, ]
  )
  
  # Print key results
  cat("\nKey Results:\n")
  cat(paste("Total cost (TCZ+GC): £", round(markov_results$total_costs_sum[1], 2), "\n"))
  cat(paste("Total QALYs (TCZ+GC):", round(markov_results$total_qalys_sum[1], 4), "\n"))
  cat(paste("Total cost (GC alone): £", round(markov_results$total_costs_sum[2], 2), "\n"))
  cat(paste("Total QALYs (GC alone):", round(markov_results$total_qalys_sum[2], 4), "\n"))
  cat(paste("ICER: £", round(markov_results$icer, 2), " per QALY\n"))
}

# Create cost effectiveness results summary table
summary_table <- data.frame(
  Scenario = character(),
  Subgroup = character(),
  Treatment_Duration_Months = numeric(),
  Total_Cost_TCZ = numeric(),
  Total_Cost_GC = numeric(),
  Inc_Cost = numeric(),
  Total_QALY_TCZ = numeric(),
  Total_QALY_GC = numeric(),
  Inc_QALY = numeric(),
  ICER = numeric(),
  stringsAsFactors = FALSE
)

for (scenario_name in names(scenario_results)) {
  res <- scenario_results[[scenario_name]]
  
  summary_table <- rbind(summary_table, data.frame(
    Scenario = scenario_name,
    Subgroup = res$scenario_info$subgroup,
    Treatment_Duration_Months = res$scenario_info$treatment_duration_months,
    Total_Cost_TCZ = round(res$markov_results$total_costs_sum[1], 2),
    Total_Cost_GC = round(res$markov_results$total_costs_sum[2], 2),
    Inc_Cost = round(res$markov_results$inc_cost, 2),
    Total_QALY_TCZ = round(res$markov_results$total_qalys_sum[1], 4),
    Total_QALY_GC = round(res$markov_results$total_qalys_sum[2], 4),
    Inc_QALY = round(res$markov_results$inc_qalys, 4),
    ICER = round(res$markov_results$icer, 2),
    stringsAsFactors = FALSE
  ))
}

# Print summary table
cat("\n========================================\n")
cat("SUMMARY OF ALL SCENARIOS\n")
cat("========================================\n\n")
print(summary_table, row.names = FALSE)

