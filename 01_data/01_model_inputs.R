#############################
# Define Model Parameters #
############################

# Intervention arms to compare
n_arms <- 2
arms <- c("TCZ_plus_GC", "GC_alone")

# Time horizon and cycle length
cycle_length <- 1 
time_horizon <- 30 * 48  # 30-year lifetime horizon on a weekly basis
n_cycles <- time_horizon / cycle_length

# Apply discount rates to costs and utilities
disc_rate <- 0.035 

apply_discount <- function(x, r, t) {
  x / ((1 + r)^t)
}

# Transition Probabilities
# Probability of First Flare (0-26 weeks)
p_first_flare_TCZ_0_26 <- 0.421
  p_first_flare_TCZ_0_26_se <- 0.016
p_first_flare_GC_0_26 <- 0.832
  p_first_flare_GC_0_26_se <- 0.015

# Probability of First Flare (27 weeks to end of GC treatment)
p_first_flare_TCZ_27_end <- 0.213
  p_first_flare_TCZ_27_end_se <- 0.006
p_first_flare_GC_27_end <- 0.987
  p_first_flare_GC_27_end_se <- 0.052

# Probability of Relapse/Flare from Remission (Escape)
p_escape_rate_TCZ <- -1.056  # log scale
  p_escape_se_TCZ <- 0.354
  p_escape_days_TCZ <- 228.39
p_escape_rate_GC <- -0.300  # log scale
  p_escape_se_GC <- 0.224
  p_escape_days_GC <- 224.77

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

# Utilities
# Baseline utilities for health states (from GiACTA trial)
utility_remission <- 0.77130
  utility_remission_alpha <- 3060.148
  utility_remission_beta <- 907.372
utility_relapse <- 0.64200
  utility_relapse_alpha <- 245.783
  utility_relapse_beta <- 137.057
  
# Utility decrements 
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
  

