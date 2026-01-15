#######################################
# RUN BASELINE ANALYSIS
#######################################

baseline_results <- run_markov_model(use_baseline = TRUE)

# Print baseline results
cat("\n=== BASELINE RESULTS ===\n")
cat(paste("Total cost (Tocilizumab + GC): £", round(baseline_results$total_costs_sum[1], 2), "\n"))
cat(paste("Total QALYs (Tocilizumab + GC):", round(baseline_results$total_qalys_sum[1], 4), "\n"))
cat(paste("Total cost (GC alone): £", round(baseline_results$total_costs_sum[2], 2), "\n"))
cat(paste("Total QALYs (GC alone):", round(baseline_results$total_qalys_sum[2], 4), "\n"))
cat(paste("Incremental cost: £", round(baseline_results$inc_cost, 2), "\n"))
cat(paste("Incremental QALYs:", round(baseline_results$inc_qalys, 4), "\n"))
cat(paste("ICER: £", round(baseline_results$icer, 2), " per QALY\n"))

# Generate and print adverse event summaries
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


#######################################
# RUN PSA ANALYSIS
#######################################

n_sim <- 10
set.seed(12345)

# Generate PSA parameters
psa_params_all <- generate_psa_params(n_sim = n_sim, seed = 12345)

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

# Run PSA simulations
for (sim in 1:n_sim) {
  
  # Extract parameters for this simulation
  sim_params <- list(
    scalar_params = as.list(psa_params_all$scalar_params[sim, ]),
    vector_params = lapply(psa_params_all$vector_params, function(x) x[sim, ])
  )
  
  # Run model with these parameters
  results <- run_markov_model(params = sim_params, use_baseline = FALSE)
  
  # Store results
  psa_total_costs[sim, ] <- results$total_costs_sum
  psa_total_qalys[sim, ] <- results$total_qalys_sum
  psa_icer[sim] <- results$icer
  
  # Store adverse event counts
  ae_summary <- summarize_adverse_events(results)
  psa_incident_events[sim, , 1] <- ae_summary$total_incident_events$TCZ_plus_GC
  psa_incident_events[sim, , 2] <- ae_summary$total_incident_events$GC_alone
  
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
plot(inc_qalys_psa, inc_costs_psa, 
     xlab = "Incremental QALYs", 
     ylab = "Incremental Costs (£)",
     main = "Cost-Effectiveness Plane",
     pch = 20, col = rgb(0, 0, 1, 0.3))
abline(h = 0, v = 0, lty = 2)

wtp_thresholds <- c(20000, 30000)
for (wtp in wtp_thresholds) {
  abline(a = 0, b = wtp, col = "red", lty = 2)
}
legend("topright", legend = c("£20,000/QALY", "£30,000/QALY"), 
       col = "red", lty = 2)

# CEAC
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

# Save results
psa_output <- data.frame(
  simulation = 1:n_sim,
  cost_tcz = psa_total_costs[, 1],
  qaly_tcz = psa_total_qalys[, 1],
  cost_control = psa_total_costs[, 2],
  qaly_control = psa_total_qalys[, 2],
  inc_cost = inc_costs_psa,
  inc_qaly = inc_qalys_psa,
  icer = psa_icer
)
write.csv(psa_output, "psa_results.csv", row.names = FALSE)

