##################################
# Markov model for diabetes #####
##################################

markov_dm <- function(l.params) {
  with(l.params, {
     
     # Define cost and QALY matrices
     ################################ 
     # cost of staying in state
     state_c_matrix_monthly <-
       matrix(c(0, c_dm/12,
                0, c_dm/12,
                0, c_dm/12,
                0, c_dm/12),
              byrow = TRUE,
              nrow = n.steroid.taper.schedule,
              dimnames = list(steroid.taper.schedule,
                              state.names))
     
     # QALY when staying in state
     state_q_matrix_monthly <-
       matrix(c(utility_remission/12, q_dm/12,
                utility_remission/12, q_dm/12,
                (utility_remission-utility_decrement_relapse)/12, (q_dm-utility_decrement_relapse)/12,
                (utility_remission-utility_decrement_relapse)/12, (q_dm-utility_decrement_relapse)/12),
              byrow = TRUE,
              nrow = n.steroid.taper.schedule,
              dimnames = list(steroid.taper.schedule,
                              state.names))
     
     # Initialize state populations
     pop <- array(data = 0,
                  dim = c(n.states, n_cycles, n.steroid.taper.schedule),
                  dimnames = list(state = state.names,
                                  cycle = NULL,
                                  diagnoses = steroid.taper.schedule))
     pop["No disease", cycle = 1, ] <- 1
     
     # Initialize storage for costs and QALYs
     cycle_empty_array <- array(NA,
                                dim = c(n.steroid.taper.schedule, n_cycles),
                                dimnames = list(schedule = steroid.taper.schedule, cycle = NULL))
     # Sum costs and QALYs for each cycle at a time for each diagnostic outcome
     cycle_costs <- cycle_qalys <- cycle_empty_array
     cycle_state_costs <- cycle_state_qalys <- cycle_empty_array
     
     # Function to calculate the time-dependent transition matrix
     p_matrix <- array(0, dim = c(n.states, n.states, n.steroid.taper.schedule),
                       dimnames = list(from = state.names, to = state.names, steroid.taper.schedule))
     
     calculate_transition_matrix <- function(type, age, cycle) {
         
         p.dm.steroid.tp.fp <- p_dm * v.steroid.hr.diabetes.tp.fp[cycle] # Steroid-related risk of developing diseases for TPs and FPs
         p.dm.steroid.tn <- p.dm * v.steroid.hr.diabetes.tn[cycle] # Steroid-related risk of developing diseases for TNs
         p.dm.steroid.fn <- p.dm * (v.steroid.hr.diabetes.fn1[cycle]*v.fn.detection[1] + 
                                      v.steroid.hr.diabetes.fn2[cycle]*v.fn.detection[2] +
                                      v.steroid.hr.diabetes.fn3[cycle]*v.fn.detection[3]) # Steroid-related excessive mortality for FNs
         
         # Convert annual probabilities into monthly probabilities
         annual_probs <- c("p.mortality", "p.mortality.steroid.tp.fp", "p.mortality.steroid.tn", "p.mortality.steroid.fn", 
                           "p.dm.steroid.tp.fp", "p.dm.steroid.tn", "p.dm.steroid.fn")
         for (prob in annual_probs) {
           monthly_prob <- 1 - exp(-(-log(1 - get(prob))) / 12)
           assign(paste0(prob, ".monthly"), monthly_prob)
         }
       
         # Matrix containing transition probabilities for TPs
         p_matrix["No disease", "Diseases", "TPs"] <- p.dm.steroid.tp.fp.monthly
         p_matrix["No disease", "Dead", "TPs"] <- p.mortality.steroid.tp.fp.monthly 
         p_matrix["Diseases", "Dead", "TPs"] <- p.mortality.steroid.tp.fp.monthly 
         
         # Matrix containing transition probabilities for FNs
         p_matrix["No disease", "Diseases", "FNs"] <- p.dm.steroid.fn.monthly
         p_matrix["No disease", "Dead", "FNs"] <- p.mortality.steroid.fn.monthly 
         p_matrix["Diseases", "Dead", "FNs"] <- p.mortality.steroid.fn.monthly 
         
         # Matrix containing transition probabilities for FPs
         p_matrix["No disease", "Diseases", "FPs"] <- p.dm.steroid.tp.fp.monthly
         p_matrix["No disease", "Dead", "FPs"] <-  p.mortality.steroid.tp.fp.monthly 
         p_matrix["Diseases", "Dead", "FPs"] <-  p.mortality.steroid.tp.fp.monthly 
         
         # Matrix containing transition probabilities for TNs
         p_matrix["No disease", "Diseases", "TNs"] <- p.dm.steroid.tn.monthly
         p_matrix["No disease", "Dead", "TNs"] <- p.mortality.steroid.tn.monthly 
         p_matrix["Diseases", "Dead", "TNs"] <- p.mortality.steroid.tn.monthly 
         
         # Add remaining probabilities to stay in the same state
         for (diag in d.names) {
           p_matrix["No disease", "No disease", diag] <- 1 - sum(p_matrix["No disease", , diag])
           p_matrix["Diseases", "Diseases", diag] <- 1 - sum(p_matrix["Diseases", , diag])
           p_matrix["Dead", "Dead", diag] <- 1
         }
         return(p_matrix)
         
       }
     
     
     # Running the model
     
     # Monthly Markov model
     for (i in 1:n.diagnoses) {
       
       age <- initial.age
       
       for (j in 2:n.monthly.cycles) {
         
         p_matrix_monthly <- calculate_transition_matrix("monthly", age, j-1)
         
         pop[, cycle = j, diagnoses = i] <-
           pop[, cycle = j - 1, diagnoses = i] %*% p_matrix_monthly[, , diagnoses = i]
         
         trans["Diseases", cycle = j, diagnoses = i] <-
           pop["No disease", cycle = j - 1, diagnoses = i] * p_matrix_monthly["No disease", "Diseases", diagnoses = i] 
         
         trans["Dead", cycle = j, diagnoses = i] <-
           pop["Diseases", cycle = j - 1, diagnoses = i] * p_matrix_monthly["Diseases", "Dead", diagnoses = i] 
         
         age <- ifelse(j<14, age, (age+1))
         
       }
     }
     
     for (i in 1:n.diagnoses) {
       
       for (j in 2:((n.monthly.cycles+1)/2)) {
         
         cycle_state_costs[i, j] <-
           (state_c_matrix_monthly[diagnoses = i, ] %*% pop[, cycle=j, diagnoses = i])
         
         cycle_costs[i, j] <- cycle_state_costs[i, j] 
         
         cycle_state_qalys[i, j] <- 
           (state_q_matrix_monthly[diagnoses = i, ] %*% pop[, cycle=j, diagnoses = i])
           
         cycle_qalys[i, j] <- cycle_state_qalys[i, j]
             
       }
     }
     
     for (i in 1:n.diagnoses) {
       
       for (j in 14:n.monthly.cycles) {
         
         cycle_state_costs[i, j] <-
           (state_c_matrix_monthly[diagnoses = i, ] %*% pop[, cycle=j, diagnoses = i])* 1/(1 + disc.rate)
         
         cycle_costs[i, j] <- cycle_state_costs[i, j] 
         
         cycle_qalys[i, j] <-
           state_q_matrix_monthly[diagnoses = i, ] %*% pop[, cycle=j, diagnoses = i]* 1/(1 + disc.rate)
       }
     }
     
    # Yearly Markov model
    for (i in 1:n.diagnoses) {
      
      age <- initial.age+2
      
      for (j in (n.monthly.cycles+1):(n.monthly.cycles+n.yearly.cycles)) {
        
        p_matrix_yearly <- calculate_transition_matrix("yearly", age, j-1)
        
        pop[, cycle = j, diagnoses = i] <-
          pop[, cycle = j - 1, diagnoses = i] %*% p_matrix_yearly[, , diagnoses = i]
        
        trans["Diseases", cycle = j, diagnoses = i] <-
          (pop["No disease", cycle = j - 1, diagnoses = i] * p_matrix_yearly["No disease", "Diseases", diagnoses = i]) 
        
        trans["Dead", cycle = j, diagnoses = i] <-
          pop["Diseases", cycle = j - 1, diagnoses = i] * p_matrix_yearly["Diseases", "Dead", diagnoses = i] 
        
        cycle_costs[diagnoses = i, cycle = j] <-
          (state_c_matrix_yearly[diagnoses = i, ] %*% pop[, cycle = j, diagnoses = i]) * 1/(1 + disc.rate)^(j - n.monthly.cycles + 1)
        
        cycle_qalys[diagnoses = i, cycle = j] <-
          (state_q_matrix_yearly[diagnoses = i, ] %*% pop[, cycle = j, diagnoses = i]) * 1/(1 + disc.rate)^(j - n.monthly.cycles + 1)
        
        age <- age + 1
      }
      
      total_diseases <- colSums(trans["Diseases", , ])
      total_dead <- colSums(trans["Dead", , ])
      total_costs <- rowSums(cycle_costs, na.rm = TRUE)
      total_qalys <- rowSums(cycle_qalys, na.rm = TRUE)
      
    }
    
    model_outputs <- list("total_diseases" = total_diseases, "total_dead" = total_dead, "total_costs" = total_costs, "total_qalys" = total_qalys)
    return(model_outputs)
      
})
}

