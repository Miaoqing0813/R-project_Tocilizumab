#############################
# Define Model Parameters #
############################

#### load packages ####
library(dplyr)
library(knitr)
library(ggplot2)
library(data.table)
library (here)
library(flexsurv)
library(SurvdigitizeR)
library(IPDfromKM)
library(survival)
library(mixcure)
library(httpgd)
library(kableExtra)
library (future)
library(doRNG)
library(doParallel)
library(darthtools)
library(survminer)
library(tm)
library(ggsurvfit)
library(zoo)
library(parallel)
library(openxlsx)

# Intervention arms to compare
n_arms <- 2
arms <- c(1, 2)   # 1 = treat (TCZ_plus_GC), 2 = control (GC_alone)

# Apply discount rates to costs and utilities
disc_rate <- 0.035 

apply_discount <- function(x, r, t) {
  x / ((1 + r)^(t - 1))
}

# Load input parameter Excel table
init.parameters <- read.csv(file = "01_data/Raw data/Init_params.csv", stringsAsFactors = F)
df.params.mean <- data.table(matrix(init.parameters$Mean, nrow = 1))
df.params.se <- data.table(matrix(init.parameters$SE, nrow = 1))
names(df.params.mean)<-init.parameters$Name
names(df.params.se)<-paste0(init.parameters$Name, "_se")

list2env(as.list(df.params.mean), .GlobalEnv)
list2env(as.list(df.params.se), .GlobalEnv)

# Mortality
initial_age = 71
allcause.mort  <- read.csv(file = "01_data/Raw data/All_cause_mortality.csv", stringsAsFactors = F)
v.allcause.mort  <- as.matrix(dplyr::select(allcause.mort,Total))

get_monthly_mort <- function(ages, mort_table) {
  annual_probs <- mort_table[pmin(round(ages), 100) + 1]
  monthly_probs <- 1 - exp(-(-log(1 - annual_probs)) / 12)
  return(monthly_probs)
}

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
  
# Utility Age Adjustment for Population Norms: assuming 50% split between men and women
population_utilities <- data.frame(
  age_group = c("70-74", "75-79", "80-84", "85-89", "90+"),
  age_lower = c(70, 75, 80, 85, 90),
  age_upper = c(74, 79, 84, 89, 100),
  female_utility = c(0.784, 0.730, 0.710, 0.666, 0.666),
  male_utility = c(0.801, 0.788, 0.767, 0.727, 0.656)
)
  
# Calculate average utility across age groups (50% male, 50% female)
population_utilities$avg_utility <- (population_utilities$female_utility + population_utilities$male_utility) / 2

# Function to get population utility for a given age
get_population_utility <- function(age) {
  # Find the age group that the age falls into
  idx <- which(age >= population_utilities$age_lower & age <= population_utilities$age_upper)
  
  # If age falls within a range, return that utility
  if (length(idx) > 0) {
    return(population_utilities$avg_utility[idx[1]])
  }

  # If age is above the highest range, use the last group
  return(population_utilities$avg_utility[nrow(population_utilities)])
}

# Function to calculate age-adjusted utility
# trial_age: mean age in trial population, using 70 rather than 69 as population norm utilities increase at age 70
# target_age: mean age in target population for PSA
# trial_utility: baseline utility measured from trial population
adjust_utility_for_age <- function(trial_utility, trial_age = 70, target_age = 71) {
  
  # Get population utilities for trial age and target age
  trial_utility_norm <- get_population_utility(trial_age)
  target_utility_norm <- get_population_utility(target_age)
  
  adjusted_utility <- (trial_utility/trial_utility_norm)* target_utility_norm
  
  return(adjusted_utility)
}

# Probability of Relapse/Flare from subsequent remission for different subgroups (following ERG's scenario analysis 6.3.1))
# Newly diagnosed patients
annual_rate_subseq_flare_newlydiag_control  <- 0.400   # Annual risk estimate
annual_rate_subseq_flare_newlydiag_control_se <- 0.319   # Standard error of the annual risk, calculated from IQR 0.21-0.64
  # Convert them into monthly probability and standard error: p_monthly = 1 - exp(-rate_annual / 12)
  probs_subseq_flare_newlydiag_control <- 1 - exp(-annual_rate_subseq_flare_newlydiag_control / 12)
  probs_subseq_flare_newlydiag_control_se <- (1 / 12) * exp(-annual_rate_subseq_flare_newlydiag_control / 12) * annual_rate_subseq_flare_newlydiag_control_se
  
annual_rate_subseq_flare_newlydiag_treat  <- 0.400 * (exp(-1.68821)/exp(-0.91629))   # Annual risk estimate
annual_rate_subseq_flare_newlydiag_treat_se <- 0.319 * (exp(-1.68821)/exp(-0.91629))  # Standard error of the annual risk, calculated from IQR 0.21-0.64
  # Convert them into monthly probability and standard error: p_monthly = 1 - exp(-rate_annual / 12)
  probs_subseq_flare_newlydiag_treat <- 1 - exp(-annual_rate_subseq_flare_newlydiag_treat / 12)
  probs_subseq_flare_newlydiag_treat_se <- (1 / 12) * exp(-annual_rate_subseq_flare_newlydiag_treat / 12) * annual_rate_subseq_flare_newlydiag_treat_se
  
# Relapsed patients
annual_rate_subseq_flare_relapsed_control  <- 0.400 * (exp(-0.37136)/exp(-0.91629))   # Annual risk estimate
annual_rate_subseq_flare_relapsed_control_se <- 0.319 * (exp(-0.37136)/exp(-0.91629))   # Standard error of the annual risk, calculated from IQR 0.21-0.64
  # Convert them into monthly probability and standard error: p_monthly = 1 - exp(-rate_annual / 12)
  probs_subseq_flare_relapsed_control <- 1 - exp(-annual_rate_subseq_flare_relapsed_control / 12)
  probs_subseq_flare_relapsed_control_se <- (1 / 12) * exp(-annual_rate_subseq_flare_relapsed_control / 12) * annual_rate_subseq_flare_relapsed_control_se
  
annual_rate_subseq_flare_relapsed_treat  <- 0.400 * (exp(-1.68821)/exp(-0.91629)) * (exp(-0.37136)/exp(-0.91629))   # Annual risk estimate
annual_rate_subseq_flare_relapsed_treat_se <- 0.319 * (exp(-1.68821)/exp(-0.91629)) * (exp(-0.37136)/exp(-0.91629))  # Standard error of the annual risk, calculated from IQR 0.21-0.64
  # Convert them into monthly probability and standard error: p_monthly = 1 - exp(-rate_annual / 12)
  probs_subseq_flare_relapsed_treat <- 1 - exp(-annual_rate_subseq_flare_relapsed_treat / 12)
  probs_subseq_flare_relapsed_treat_se <- (1 / 12) * exp(-annual_rate_subseq_flare_relapsed_treat / 12) * annual_rate_subseq_flare_relapsed_treat_se
  
# ITT (Intent-to-Treat) population
annual_rate_subseq_flare_itt_control  <- 0.400 * exp(-1.36928)/exp(-1.68821)   # Annual risk estimate
annual_rate_subseq_flare_itt_control_se <- 0.319 * exp(-1.36928)/exp(-1.68821)   # Standard error of the annual risk, calculated from IQR 0.21-0.64
  # Convert them into monthly probability and standard error: p_monthly = 1 - exp(-rate_annual / 12)
  probs_subseq_flare_itt_control <- 1 - exp(-annual_rate_subseq_flare_itt_control / 12)
  probs_subseq_flare_itt_control_se <- (1 / 12) * exp(-annual_rate_subseq_flare_itt_control / 12) * annual_rate_subseq_flare_itt_control_se
  
annual_rate_subseq_flare_itt_treat  <- 0.400 * (exp(-1.68821)/exp(-0.91629)) * exp(-1.36928)/exp(-1.68821)   # Annual risk estimate
annual_rate_subseq_flare_itt_treat_se <- 0.319 * (exp(-1.68821)/exp(-0.91629)) * exp(-1.36928)/exp(-1.68821)  # Standard error of the annual risk, calculated from IQR 0.21-0.64
  # Convert them into monthly probability and standard error: p_monthly = 1 - exp(-rate_annual / 12)
  probs_subseq_flare_itt_treat <- 1 - exp(-annual_rate_subseq_flare_itt_treat / 12)
  probs_subseq_flare_itt_treat_se <- (1 / 12) * exp(-annual_rate_subseq_flare_itt_treat / 12) * annual_rate_subseq_flare_itt_treat_se
  