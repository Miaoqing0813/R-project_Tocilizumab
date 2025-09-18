######################################
# Markov model for vision loss #####
######################################

markov_vl <- function(l.params) {
  with(l.params, {
    
    # Define cost and QALY matrices
    ################################
    # cost of staying in state: monthly cost for years 1 and 2 
    # Cost of vision loss is only applied for TPs and FNs as it is a GCA-related complication
    state_c_matrix_monthly <-
      matrix(c(0, c.vl.y2on/12, 0,
               0, c.vl.y2on/12, 0,
               0, 0, 0,
               0, 0, 0),
             byrow = TRUE,
             nrow = n.diagnoses,
             dimnames = list(d.names,
                             state.names))
    
    # cost of staying in state: annual costs from year 3
    state_c_matrix_yearly <-
      matrix(c(0, c.vl.y2on, 0,
               0, c.vl.y2on, 0,
               0, 0, 0, 
               0, 0, 0),
             byrow = TRUE,
             nrow = n.diagnoses,
             dimnames = list(d.names,
                             state.names))
    
    # cost of moving to a state: 
    # to account for higher cost in year 1; all difference was added on month 1
    trans_c_matrix <-
      matrix(c(0, (c.vl.y1-c.vl.y2on), 0,
               0, (c.vl.y1-c.vl.y2on), 0,
               0, 0, 0,
               0, 0, 0),
              byrow = TRUE,
              nrow = n.diagnoses,
              dimnames = list(d.names,
                              state.names))
    
    # QALY when staying in state: monthly QALY for years 1 and 2
    # QALY of vision loss is only applied for TPs and FNs as it is a GCA-related complication
    state_q_matrix_monthly <-
      matrix(c(q.baseline/12, q.vl/12, 0,
               q.baseline/12, q.vl/12, 0,
               q.baseline/12, q.vl/12, 0, 
               q.baseline/12, q.vl/12, 0),
             byrow = TRUE,
             nrow = n.diagnoses,
             dimnames = list(d.names,
                             state.names))
    
    # QALY when staying in state: annual utility from year 3
    state_q_matrix_yearly <-
      matrix(c(q.baseline, q.vl, 0,
               q.baseline, q.vl, 0,
               q.baseline, q.vl, 0, 
               q.baseline, q.vl, 0),
             byrow = TRUE,
             nrow = n.diagnoses,
             dimnames = list(d.names,
                             state.names))

    
    # Initialize state populations
    pop <- array(data = 0,
                 dim = c(n.states, (n.monthly.cycles+n.yearly.cycles), n.diagnoses),
                 dimnames = list(state = state.names,
                                 cycle = NULL,
                                 diagnoses = d.names))
    pop["No disease", cycle = 1, ] <- 1
    
    # People in transition
    trans <- array(data = 0,
                   dim = c(n.states, (n.monthly.cycles+n.yearly.cycles), n.diagnoses),
                   dimnames = list(state = state.names,
                                   cycle = NULL,
                                   diagnoses = d.names))
    
   
    # Initialize storage for costs and QALYs
    cycle_empty_array <- array(NA,
                               dim = c(n.diagnoses, (n.monthly.cycles+n.yearly.cycles)),
                               dimnames = list(diagnoses = d.names, cycle = NULL))
    # Sum costs and QALYs for each cycle at a time for each diagnostic outcome
    cycle_costs <- cycle_qalys <- cycle_empty_array
    cycle_state_costs <- cycle_trans_costs <- cycle_empty_array
    
    # Function to calculate the time-dependent transition matrix
    p_matrix <- array(0, dim = c(n.states, n.states, n.diagnoses),
                      dimnames = list(from = state.names, to = state.names, d.names))
    
    calculate_transition_matrix <- function(type, age, cycle) {
      if (type == "monthly") {
        
        p.mortality <- v.allcause.mort[age+1] # Age-dependent mortality
        p.mortality.steroid.tp.fp <- p.mortality * v.steroid.hr.mortality.tp.fp[cycle] # Steroid-related excessive mortality for TPs and FPs
        p.mortality.steroid.tn <- p.mortality * v.steroid.hr.mortality.tn[cycle] # Steroid-related excessive mortality for TNs 
        p.mortality.steroid.fn <- p.mortality * (v.steroid.hr.mortality.fn1[cycle]*v.fn.detection[1] + 
                                                   v.steroid.hr.mortality.fn2[cycle]*v.fn.detection[2] +
                                                   v.steroid.hr.mortality.fn3[cycle]*v.fn.detection[3]) # Steroid-related excessive mortality for FNs
        
        v.p.vl.fn <- c(p.vl.fn, # month 1: everyone is FN
                       p.vl.fn*(1-v.fn.detection[1])+p.vl.tp*v.fn.detection[1], # Month2: 25% are detected
                       p.vl.fn*(1-v.fn.detection[1]-v.fn.detection[2])+p.vl.tp*(v.fn.detection[1]+v.fn.detection[2]), # Month3: another 50% are detected
                       rep(p.vl.tp,13), # From month 4 onwards, everyone has become TPs
                       rep(0, 9)) # Assume the incidence of vision loss is only developed in year 1 after steroid use
        p.vl.fn <- v.p.vl.fn[cycle]
        
        # Assume the incidence of vision loss is only developed in year 1
        v.p.vl.tp <- c(rep(p.vl.tp, 13), rep(0, 12))
        p.vl.tp <- v.p.vl.tp[cycle]
        
        # Convert annual probabilities into monthly probabilities
        annual_probs <- c("p.mortality", "p.mortality.steroid.tp.fp", "p.mortality.steroid.tn", "p.mortality.steroid.fn", "p.vl.tp", "p.vl.fn")
        for (prob in annual_probs) {
          monthly_prob <- 1 - exp(-(-log(1 - get(prob))) / 12)
          assign(paste0(prob, ".monthly"), monthly_prob)
        }
        
        # Matrix containing transition probabilities for TPs
        p_matrix["No disease", "Diseases", "TPs"] <- p.vl.tp.monthly
        p_matrix["No disease", "Dead", "TPs"] <- p.mortality.steroid.tp.fp.monthly 
        p_matrix["Diseases", "Dead", "TPs"] <- p.mortality.steroid.tp.fp.monthly 
        
        # Matrix containing transition probabilities for FNs
        p_matrix["No disease", "Diseases", "FNs"] <- p.vl.fn.monthly
        p_matrix["No disease", "Dead", "FNs"] <- p.mortality.steroid.fn.monthly 
        p_matrix["Diseases", "Dead", "FNs"] <- p.mortality.steroid.fn.monthly 
        
        # Matrix containing transition probabilities for FPs
        p_matrix["No disease", "Diseases", "FPs"] <- 0
        p_matrix["No disease", "Dead", "FPs"] <- p.mortality.steroid.tp.fp.monthly 
        p_matrix["Diseases", "Dead", "FPs"] <- p.mortality.steroid.tp.fp.monthly 
        
        # Matrix containing transition probabilities for TNs
        p_matrix["No disease", "Diseases", "TNs"] <- 0
        p_matrix["No disease", "Dead", "TNs"] <- p.mortality.steroid.tn.monthly
        p_matrix["Diseases", "Dead", "TNs"] <- p.mortality.steroid.tn.monthly
        
        # Add remaining probabilities to stay in the same state
        for (diag in d.names) {
          p_matrix["No disease", "No disease", diag] <- 1 - sum(p_matrix["No disease", , diag])
          p_matrix["Diseases", "Diseases", diag] <- 1 - sum(p_matrix["Diseases", , diag])
          p_matrix["Dead", "Dead", diag] <- 1
        }
        return(p_matrix)
        
      } else if (type == "yearly") {
        
        p.mortality <- v.allcause.mort[age+1] # Age-dependent mortality
        
        # Adjust for the length of steroid effects: 
        v.p.mortality.steroid.tp.fp <- c(rep(p.mortality*v.steroid.hr.mortality.tp.fp[n.monthly.cycles], n.monthly.cycles+steroid.effect.length-1), # Increased mortality due to steroid use
                                         rep(p.mortality, n.yearly.cycles-steroid.effect.length+1)) # Mortality in the general population for the rest of time horizon
        p.mortality.steroid.tp.fp<- v.p.mortality.steroid.tp.fp[cycle]
        v.p.mortality.steroid.tn <- c(rep(p.mortality*v.steroid.hr.mortality.tn[n.monthly.cycles], n.monthly.cycles+steroid.effect.length-2), # Increased mortality due to steroid use: minus one year for TNs
                                      rep(p.mortality, n.yearly.cycles-steroid.effect.length+2)) # Mortality in the general population for the rest of time horizon
        p.mortality.steroid.tn<- v.p.mortality.steroid.tn[cycle]
        
        # Matrix containing transition probabilities for TPs
        p_matrix["No disease", "Diseases", "TPs"] <- 0
        p_matrix["No disease", "Dead", "TPs"] <- p.mortality.steroid.tp.fp 
        p_matrix["Diseases", "Dead", "TPs"] <- p.mortality.steroid.tp.fp 
        
        # Matrix containing transition probabilities for FNs: same as TPs 
        p_matrix["No disease", "Diseases", "FNs"] <- 0
        p_matrix["No disease", "Dead", "FNs"] <- p.mortality.steroid.tp.fp 
        p_matrix["Diseases", "Dead", "FNs"] <- p.mortality.steroid.tp.fp 
        
        # Matrix containing transition probabilities for FPs
        p_matrix["No disease", "Diseases", "FPs"] <- 0
        p_matrix["No disease", "Dead", "FPs"] <- p.mortality.steroid.tp.fp 
        p_matrix["Diseases", "Dead", "FPs"] <- p.mortality.steroid.tp.fp 
        
        # Matrix containing transition probabilities for TNs
        p_matrix["No disease", "Diseases", "TNs"] <- 0
        p_matrix["No disease", "Dead", "TNs"] <- p.mortality.steroid.tn 
        p_matrix["Diseases", "Dead", "TNs"] <- p.mortality.steroid.tn 
        
        # Add remaining probabilities to stay in the same state
        for (diag in d.names) {
          p_matrix["No disease", "No disease", diag] <- 1 - sum(p_matrix["No disease", , diag])
          p_matrix["Diseases", "Diseases", diag] <- 1 - sum(p_matrix["Diseases", , diag])
          p_matrix["Dead", "Dead", diag] <- 1
        }
        return(p_matrix)
        
      }
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
        
        cycle_trans_costs[i, j] <- trans["Diseases", cycle = j, diagnoses = i] * trans_c_matrix[diagnoses = i,"Diseases"]
        
        cycle_costs[i, j] <- cycle_state_costs[i, j] + cycle_trans_costs[i, j]
        
        cycle_qalys[i, j] <-
          state_q_matrix_monthly[diagnoses = i, ] %*% pop[, cycle=j, diagnoses = i]
      }
    }
    
    for (i in 1:n.diagnoses) {
      
      for (j in 14:n.monthly.cycles) {
        
        cycle_state_costs[i, j] <-
          (state_c_matrix_monthly[diagnoses = i, ] %*% pop[, cycle=j, diagnoses = i])* 1/(1 + disc.rate)
        
        cycle_trans_costs[i, j] <- trans["Diseases", cycle = j, diagnoses = i] * trans_c_matrix[diagnoses = i,"Diseases"]*1/(1 + disc.rate)
        
        cycle_costs[i, j] <- cycle_state_costs[i, j] + cycle_trans_costs[i, j]
        
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

        cycle_state_costs[diagnoses = i, cycle = j] <-
          (state_c_matrix_yearly[diagnoses = i, ] %*% pop[, cycle = j, diagnoses = i]) * 1/(1 + disc.rate)^(j - n.monthly.cycles + 1)
        
        cycle_trans_costs[diagnoses = i, cycle = j] <-
          trans["Diseases", cycle = j, diagnoses = i] * trans_c_matrix[diagnoses = i,"Diseases"] * 1/(1 + disc.rate)^(j - n.monthly.cycles + 1) 
        
        cycle_costs[i, j] <- cycle_state_costs[i,j] + cycle_trans_costs[i, j]
        
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



  



