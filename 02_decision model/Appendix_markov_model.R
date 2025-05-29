# --------------------- 4. Summaries ------------------------------------

# Run Markov model
for (cycle in 2:n_cycles) {
  for (arm in arms) {
    prev <- pop[, cycle - 1, arm]
    TPM <- ifelse(arm == "TCZ_plus_GC", tpm_TCZ, tpm_GC)
    pop[, cycle, arm] <- as.numeric(prev %*% TPM)
  }
}

# Assign utility and cost vectors
utility_vector <- c(utility_remission_new_onset,
                    utility_flare_new_onset,
                    utility_remission_relapse,
                    0)

cost_vector_TCZ <- c(cost_TCZ_PAS,
                     cost_diabetes_per_week + cost_fracture,
                     cost_TCZ_PAS,
                     0)

cost_vector_GC <- c(cost_GC_per_week,
                    cost_diabetes_per_week + cost_fracture,
                    cost_GC_per_week,
                    0)

# Cycle-by-cycle accumulation
for (cycle in 1:n_cycles) {
  for (arm in arms) {
    cost_vector <- if (arm == "TCZ_plus_GC") cost_vector_TCZ else cost_vector_GC
    state_dist <- pop[, cycle, arm]
    
    # Apply discounting
    discounted_costs <- apply_discount(cost_vector, disc_rate, t = cycle)
    discounted_utils <- apply_discount(utility_vector, disc_rate, t = cycle)
    
    cycle_costs[arm, cycle] <- sum(state_dist * discounted_costs)
    cycle_qalys[arm, cycle] <- sum(state_dist * discounted_utils)
  }
}

# Summarise total outcomes
total_costs <- colSums(cycle_costs)
total_qalys <- colSums(cycle_qalys)


total_cost_TCZ_plus <- sum(costs_TCZ_plus)
total_qalys_TCZ_plus<- sum(qalys_TCZ_plus)

total_cost_GC_alone <- sum(costs_GC_alone)
total_qalys_GC_alone<- sum(qalys_GC_alone)

# Incremental
inc_cost  <- total_cost_TCZ_plus - total_cost_GC_alone
inc_qalys <- total_qalys_TCZ_plus - total_qalys_GC_alone

icer <- inc_cost / inc_qalys

# Print results
cat("==== Tocilizumab plus GC vs. GC Alone ====\n")
cat("Total Costs (TCZ + GC):   ", round(total_cost_TCZ_plus,2), "\n")
cat("Total QALYs (TCZ + GC):  ", round(total_qalys_TCZ_plus,2), "\n\n")

cat("Total Costs (GC Alone):  ", round(total_cost_GC_alone,2), "\n")
cat("Total QALYs (GC Alone): ",  round(total_qalys_GC_alone,2), "\n\n")

cat("Incremental Cost:        ", round(inc_cost,2), "\n")
cat("Incremental QALYs:       ", round(inc_qalys,2), "\n")
cat("ICER:                    ", round(icer,2), "\n")

