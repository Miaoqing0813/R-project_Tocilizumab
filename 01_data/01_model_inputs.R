#############################
# Define Model Parameters #
############################

# Intervention arms to compare
n_arms <- 2
arms <- c(1, 0)   # 1 = treat (TCZ_plus_GC), 0 = control (GC_alone)

# Time horizon and cycle length
cycle_length <- 1 
n_cycles <- 10 * 12  # 10-year lifetime horizon on a monthly basis
time_horizon_weeks <- 10*52 # 10-year lifetime horizon on a weekly basis

# Patient's age 
initial_age <- 71 # starting age
v.age <- initial_age + ((1:n_cycles - 1) * cycle_length / 12) # get patients' age in each cycle

# Steroid effects duration
l.steroid.effect.durations <- c("5y" = 5*12, 
                                "10y" = 10*12, 
                                "15y" = 15*12, 
                                "20y" = 20*12, 
                                "30y" = 30*12)

steroid_effect_duration_months <- l.steroid.effect.durations["5y"] # base case for 5 years

# Apply discount rates to costs and utilities
disc_rate <- 0.035 

apply_discount <- function(x, r, t) {
  x / ((1 + r)^(t - 1))
}

# Different steroid tapering schedules
steroid.taper.schedule <- c("Control", "Treatment", "Control_relapse", "Tretment_relapse") 
n.steroid.taper.schedule <- length(steroid.taper.schedule)

# Health states of Markov models for steroid-related adverse events
state.names <- c("No disease", "Diseases") # a 2-state Markov model
n.states <- length(state.names)

# Load input parameter Excel table
init.parameters <- read.csv(file = "01_data/Raw data/Init_params.csv", stringsAsFactors = F)
df.params.mean <- data.table(matrix(init.parameters$Mean, nrow = 1))
names(df.params.mean)<-init.parameters$Name
list2env(as.list(df.params.mean), .GlobalEnv)

# Mortality
allcause.mort  <- read.csv(file = "01_data/Raw data/All_cause_mortality.csv", stringsAsFactors = F)
v.allcause.mort  <- as.matrix(dplyr::select(allcause.mort,Weighted))

get_monthly_mort <- function(ages, mort_table) {
  annual_probs <- mort_table[pmin(round(ages), 100) + 1]
  monthly_probs <- 1 - exp(-(-log(1 - annual_probs)) / 12)
  return(monthly_probs)
}
v.allcause.mort.monthly <- get_monthly_mort(v.age, v.allcause.mort)

# Hazard ratios for steroid-related adverse events
steroid.hr  <- read.csv(file = "01_data/Raw data/Steroid HR.csv", stringsAsFactors = F)
  # Define the variable groups and suffixes
  conditions <- c("Mortality", "Diabetes", "CVD", "Infection", "Fractures")
  group <- c("Control", "Treat", "Control.Relapse", "Treat.Relapse")
  suffixes <- c("", ".LCI", ".UCI")  # For mean, LCI, and UCI

  # Loop through each combination of conditions, treatment groups, and suffix
  for (var in conditions) {
    for (path in group) {
      for (suffix in suffixes) {
        # Construct the column name dynamically
        column_name <- paste0(var, ".", path, suffix)
        
        # Construct the variable name dynamically
        var_name <- paste0("v.steroid.hr.", tolower(var), ".", tolower(path), tolower(suffix))
        
        # Assign the selected column to a matrix variable dynamically
        assign(var_name, as.matrix(dplyr::select(steroid.hr, !!sym(column_name))))
      }
    }
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
      stroke     = p_stroke,
      mi         = p_mi,
      hf         = p_hf,
      dm         = p_dm,
      infection  = p_infection,
      fractures  = p_fractures,
      mortality  = v.allcause.mort.monthly[cycle]
    )
    
  } else {
    if (state == "Initial remission") {
      idx <- min(cycle, max_cycle)
    } else if (state == "Relapse") {
      idx <- 1
    } else if (state == "Subsequent remission") {
      # Use cycles since entering this state (1st cycle in this state = index 2)
      idx <- min(time_in_state + 1, max_cycle)
    }
    
    # Select HRs by schedule
    if (schedule_index == "treat") {
      cvd_hr <- v.steroid.hr.cvd.treat[idx, ]
      dm_hr <- v.steroid.hr.diabetes.treat[idx, ]
      infection_hr <- v.steroid.hr.infection.treat[idx, ]
      fractures_hr <- v.steroid.hr.fractures.treat[idx, ]
      mortality_hr <- v.steroid.hr.mortality.treat[idx, ]
    } else if (schedule_index == "control") {
      cvd_hr <- v.steroid.hr.cvd.control[idx, ]
      dm_hr <- v.steroid.hr.diabetes.control[idx, ]
      infection_hr <- v.steroid.hr.infection.control[idx, ]
      fractures_hr <- v.steroid.hr.fractures.control[idx, ]
      mortality_hr <- v.steroid.hr.mortality.control[idx, ]
    } else if (schedule_index == "treat.relapse") {
      cvd_hr <- v.steroid.hr.cvd.treat.relapse[idx, ]
      dm_hr <- v.steroid.hr.diabetes.treat.relapse[idx, ]
      infection_hr <- v.steroid.hr.infection.treat.relapse[idx, ]
      fractures_hr <- v.steroid.hr.fractures.treat.relapse[idx, ]
      mortality_hr <- v.steroid.hr.mortality.treat.relapse[idx, ]
    } else if (schedule_index == "control.relapse") {
      cvd_hr <- v.steroid.hr.cvd.control.relapse[idx, ]
      dm_hr <- v.steroid.hr.diabetes.control.relapse[idx, ]
      infection_hr <- v.steroid.hr.infection.control.relapse[idx, ]
      fractures_hr <- v.steroid.hr.fractures.control.relapse[idx, ]
      mortality_hr <- v.steroid.hr.mortality.control.relapse[idx, ]
    }
    
    steroid_ae_probs <- list(
      stroke     = p_stroke    * cvd_hr,
      mi         = p_mi        * cvd_hr,
      hf         = p_hf        * cvd_hr,
      dm         = p_dm        * dm_hr,
      infection  = p_infection * infection_hr,
      fractures  = p_fractures * fractures_hr,
      mortality  = v.allcause.mort.monthly[cycle] * mortality_hr
    )
  }
  
  return(steroid_ae_probs)
}


# Probability of First Flare
source("01_data/02b_survival subgroups.R")
trans_probs_first_flare_treat_weibull # cycle-specific first flare probabilities for the TCZ arm
trans_probs_first_flare_control_weibull # cycle-specific first flare probabilities for the control arm
trans_probs_first_flare_treat_exp # cycle-specific first flare probabilities for the TCZ arm
trans_probs_first_flare_control_exp # cycle-specific first flare probabilities for the control arm

# Probability of Relapse/Flare from subsequent remission (Escape) for newly-onset patients following ERG's scenario analysis 6.3.1
probs_subseq_flare_treat <- 0.0035*(52/12)
probs_subseq_flare_control <- 0.0076*(52/12)

# Utilities
# Utilities for health states (ITT population)
utility_remission <- 0.77130
  utility_remission_alpha <- 3060.148
  utility_remission_beta <- 907.372
utility_relapse <- 0.64200
  utility_relapse_alpha <- 245.783
  utility_relapse_beta <- 137.057
utility_decrement_relapse <- utility_remission - utility_relapse
  
# Utility decrements 
utility_decrement_vl <- utility_remission - q_vl
utility_decrement_stroke_y1 <- utility_remission - q_stroke_y1
utility_decrement_stroke_y2on <- utility_remission - q_stroke_y1
utility_decrement_mi_y1 <- utility_remission - q_mi_y1
utility_decrement_mi_y2on <- utility_remission - q_mi_y1
utility_decrement_hf_y1 <- utility_remission - q_hf_y1
utility_decrement_hf_y2on <- utility_remission - q_hf_y1
utility_decrement_fracture_y1 <- utility_remission - q_fracture_y1
utility_decrement_fracture_y2on <- utility_remission - q_fracture_y1
utility_decrement_dm <- utility_remission - q_dm

utility_decrement_GCA <- 0.1293
utility_decrement_GC <- 0.03
utility_decrement_vision_loss <- 0.36734
  utility_decrement_vision_loss_alpha <- 217.927
  utility_decrement_vision_loss_beta <- 0.003
utility_decrement_minor_stroke <- 0.17882
  utility_decrement_minor_stroke_alpha <- 754.201
  utility_decrement_minor_stroke_beta <- 0.001
utility_decrement_major_stroke <- 0.49122
  utility_decrement_major_stroke_alpha <- 105.394
  utility_decrement_major_stroke_beta <- 0.005
utility_decrement_diabetes <- 0.09264
  utility_decrement_diabetes_alpha <- 177.372
  utility_decrement_diabetes_beta <- 0.001
utility_decrement_fractures_vertebral_y1 <- 0.33179
  utility_decrement_fractures_vertebral_y1_alpha <- 269.152
  utility_decrement_fractures_vertebral_y1_beta <- 0.002
utility_decrement_fractures_hip_y1 <- 0.23915
  utility_decrement_fractures_hip_y1_alpha <- 484.138
  utility_decrement_fractures_hip_y1_beta <- 0.002
utility_decrement_fractures_humerus_y1 <- 0.10772
  utility_decrement_fractures_humerus_y1_alpha <- 1478.150
  utility_decrement_fractures_humerus_y1_beta <- 0.001
utility_decrement_fractures_forearm_y1 <- 0.09264
  utility_decrement_fractures_forearm_y1_alpha <- 177.372
  utility_decrement_fractures_forearm_y1_beta <- 0.001
utility_decrement_fractures_vertebral_y2on <- 0.26177
  utility_decrement_fractures_vertebral_y2on_alpha <- 416.389
  utility_decrement_fractures_vertebral_y2on_beta <- 0.002
utility_decrement_fractures_hip_y2on <- 0.11526
  utility_decrement_fractures_hip_y2on_alpha <- 1358.198
  utility_decrement_fractures_hip_y2on_beta <- 0.001
utility_decrement_fractures_humerus_y2on <- 0.00000
  utility_decrement_fractures_humerus_y2on_alpha <- 0.000
  utility_decrement_fractures_humerus_y2on_beta <- 0.000
utility_decrement_fractures_forearm_y2on <- 0.01508
  utility_decrement_fractures_forearm_y2on_alpha <- 12864.486
  utility_decrement_fractures_forearm_y2on_beta <- 0.000
  
# Costs
# Tocilizumab (TCZ) Costs: 4*162 mg syringe
cost_TCZ <- 913.12   
cost_TCZ_PAS <- cost_TCZ * 0.75  # Applying a 25% discount for PAS (adjust as needed)

# Prednisolone (GC) Costs: 30 tablets at 5mg each
cost_GC_per_week <- 26.7  

# Costs for Adverse Events
cost_diabetes <- 48.30
  cost_diabetes_se <- 0.15
cost_fractures_vertebral <- 1152
  cost_fractures_vertebral_se <- 0.15
cost_fractures_hip <- 4222
  cost_fractures_hip_se <- 0.15
cost_fractures_forearm <- 690
  cost_fractures_forearm_se <- 0.15
cost_fractures_humerus <- 690
  cost_fractures_humerus_se <- 0.15
  
cost_vision_loss_y1 <- 97.55
  cost_vision_loss_y1_se <- 0.15
cost_vision_loss_y2on <- 93.97
  cost_vision_loss_y2on_se <- 0.15
cost_nonfatal_stroke <- 112.69
  cost_nonfatal_stroke_se <- 0.15
  
# Treatment costs in primary care
cost_primary_care <- 115.56
  cost_primary_care_se <- 0.15

# Cost of primary care referral to specialist care
cost_referral <- 568.88
  cost_referral_se <- 0.15
  
################################################################################
###### Input parameters used in HTA model #####################################
###############################################################################
 
# GCA-Related Complications at Relapse/Flare
p_visual_loss <- 0.00025
p_visual_loss_alpha <- 12
p_visual_loss_beta <- 900
p_stroke <- 0.00050
p_minor_stroke <- 0.0030
p_major_stroke <- 0.0020
p_stroke_alpha <- 60
p_stroke_beta <- 895

# Steroid-Related Adverse Events
p_fractures_intercept <- -2.6583
p_fractures_intercept_se <- 0.0849
p_fractures_slope <- 0.0746
p_fractures_intercept_se <- 0.0275
p_diabetes_intercept <- -3.9216
p_diabetes_intercept_se <- 0.1346
p_diabetes_slope <- 0.1297
p_diabetes_intercept_se <- 0.0373