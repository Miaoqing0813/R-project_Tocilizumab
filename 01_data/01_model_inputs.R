#############################
# Define Model Parameters #
############################

# Intervention arms to compare
n_arms <- 2
arms <- c(1, 0)   # 1 = treat (TCZ_plus_GC), 0 = control (GC_alone)

# Time horizon and cycle length
cycle_length <- 1 
n_cycles <- 5 * 12  # 5-year lifetime horizon on a monthly basis
time_horizon_weeks <- 5 * 52 # 5-year lifetime horizon on a weekly basis

# Patient's age 
initial_age <- 71 # starting age
v.age <- initial_age + ((1:n_cycles - 1) * cycle_length / 12) # get patients' age in each cycle

# Steroid effects duration
steroid_effect_duration_months <- 5*12

# Apply discount rates to costs and utilities
disc_rate <- 0.035 

apply_discount <- function(x, r, t) {
  x / ((1 + r)^(t - 1))
}

# Different steroid tapering schedules depending on health states and treatment arms
steroid.taper.schedule <- c("Control", "Treatment", "Control_relapse", "Tretment_relapse") 

# Load input parameter Excel table
init.parameters <- read.csv(file = "01_data/Raw data/Init_params.csv", stringsAsFactors = F)
df.params.mean <- data.table(matrix(init.parameters$Mean, nrow = 1))
df.params.se <- data.table(matrix(init.parameters$SE, nrow = 1))
names(df.params.mean)<-init.parameters$Name
names(df.params.se)<-paste0(init.parameters$Name, "_se")

list2env(as.list(df.params.mean), .GlobalEnv)
list2env(as.list(df.params.se), .GlobalEnv)

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
steroid_events <- c("stroke", "mi", "hf", "fractures", "dm", "infection")
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

# TCZ monthly costing for first 12 cycles (one year) 
tcz_months <- 12L # assume one-year TCZ treatment; notice that the monitoring costs for the steroid mono therapy arm is assumed for 12 months as well
is_on_tcz <- function(cycle) cycle <= tcz_months

# Probability of First Flare for each cycle
source("01_data/02b_survival subgroups.R")
trans_probs_first_flare_treat_weibull 
trans_probs_first_flare_control_weibull 
trans_probs_first_flare_treat_exp 
trans_probs_first_flare_control_exp 

# Probability of Relapse/Flare from subsequent remission (Escape) for newly-onset patients following ERG's scenario analysis 6.3.1
probs_subseq_flare_treat <- 0.0035*(52/12)
probs_subseq_flare_control <- 0.0076*(52/12)
  

