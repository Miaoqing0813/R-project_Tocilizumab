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

# Treatment duration 
is_on_tcz <- function(cycle) cycle <= treatment_duration_months

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
allcause.mort  <- read.csv(file = "01_data/Raw data/All_cause_mortality.csv", stringsAsFactors = F)
v.allcause.mort  <- as.matrix(dplyr::select(allcause.mort,Weighted))

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

# Probability of Relapse/Flare from subsequent remission for different subgroups (following ERG's scenario analysis 6.3.1))
# Newly diagnosed patients 
probs_subseq_flare_newlydiag_treat <- 0.0035*(52/12)
probs_subseq_flare_newlydiag_control <- 0.0076*(52/12)

# Relapsed patients
probs_subseq_flare_relapsed_treat <- 0.0061*(52/12)
probs_subseq_flare_relapsed_control <- 0.0131*(52/12)

# ITT (Intent-to-Treat) population
probs_subseq_flare_itt_treat <- 0.0049*(52/12)
probs_subseq_flare_itt_control <- 0.0105*(52/12)

