######################################################################
# Time to first flare transition probability (newly-onset population)#
######################################################################

# To digitize from published Kaplan–Meier (KM) curves
first_flare_data <- survival_digitize(img_path = here(getwd(),"01_data/Raw data/Remission to first flare_subgroups_clean.png"), 
                                      num_curves = 4,
                                      x_start = 0,
                                      x_end = 52,
                                      x_increment = 4,
                                      y_start = 0,
                                      y_end = 100,
                                      y_increment = 20,
                                      y_text_vertical =  T,
                                      censoring = F)

# Generating a graph and checking how it looks compared to the original
first_flare_graph <- ggplot(data=first_flare_data) + geom_line(aes(x=time,y=St, colour=factor(curve))) + expand_limits(x = 0, y = 0) + scale_y_continuous(breaks=seq(0,100,20)) + scale_x_continuous(breaks=seq(0,52,4))
first_flare_graph

# Keeping the Km-curves generated for Tocilizumab QW + 26-week of tapering steroids and placebo + 52-week of tapering steroids for the on-set group
first_flare_data_treat <- first_flare_data[ which(first_flare_data$curve==1),c(2,3) ] # treatment group
first_flare_data_control <- first_flare_data[ which(first_flare_data$curve==4),c(2,3) ] # control group

# Generating the vector variables which provide information about the follow-up time 
# and the number of patients at risk at every time point. Both of these vectors
# will be used for the preprocess state for generating IPD datasets.
weeks <- c(0,4,8,12,16,20,24,28,32,36,40,44,48)
n_risk_treat <- c(47,43,42,40,40,40,38,36,35,34,33,33,32)
n_risk_control <- c(23,21,21,19,18,17,16,14,14,13,13,11,10)

#Preprocess state of the KM curve information 
pre_first_flare_data_treat <- preprocess(first_flare_data_treat, trisk =weeks, nrisk=n_risk_treat,totalpts = NULL,maxy=100)
pre_first_flare_data_control <- preprocess(first_flare_data_control, trisk =weeks, nrisk=n_risk_control,totalpts = NULL,maxy=100)

# Dataset IPD generation. This dataset will be used for the time to event analysis with the flexsurv package
est_first_flare_data_treat <- getIPD(prep=pre_first_flare_data_treat,armID=1)
est_first_flare_data_control <- getIPD(prep=pre_first_flare_data_control,armID=0)
IPD_first_flare_data_treat <- est_first_flare_data_treat$IPD
IPD_first_flare_data_control <- est_first_flare_data_control$IPD

# KM-curve -- within trial period
survfit2(Surv(time = time,  event = status)  ~ 1, data = IPD_first_flare_data_treat) %>%
  ggsurvfit() +
  labs(
    x = "Weeks",
    y = "Probability of 1st flare for the Tocilizumab group"
  ) + 
  add_confidence_interval() +
  add_risktable()

survfit2(Surv(time = time,  event = status)  ~ 1, data = IPD_first_flare_data_control) %>%
  ggsurvfit() +
  labs(
    x = "Weeks",
    y = "Probability of 1st flare for the comparison group"
  ) + 
  add_confidence_interval() +
  add_risktable()

# Estimating x-year survival and median survival time
km_treat <- survfit(Surv(time, status) ~ 1, data = IPD_first_flare_data_treat) 
km_control <- survfit(Surv(time, status) ~ 1, data = IPD_first_flare_data_control)  
summary(survfit(Surv(time, status) ~ 1, data = IPD_first_flare_data_treat), times = 52)
summary(survfit(Surv(time, status) ~ 1, data = IPD_first_flare_data_control), times = 52)

# Fit multiple survival models on within-trial data
surv_first_flare_treat_within <- fit.fun(
  time = "time",
  status = "status",
  data = IPD_first_flare_data_treat,
  xlow = 0, xhigh = 52,
  ylow = 0.70, yhigh = 1,
  mods = c("exp", "weibull", "gamma", "lnorm")
)
surv_first_flare_treat_within$AIC # Exponential
surv_first_flare_treat_within$BIC # Exponential
exp_model_treat_within <- surv_first_flare_treat_within$model.objects$models[["Exponential"]] 

surv_first_flare_control_within <- fit.fun(
  time = "time", 
  status = "status", 
  data = IPD_first_flare_data_control,
  xlow = 0, xhigh = 52,
  ylow = 0.40, yhigh = 1,
  mods = c("exp", "weibull", "gamma", "lnorm")
) 
surv_first_flare_control_within$AIC # Exponential
surv_first_flare_control_within$BIC # Exponential
exp_model_control_within <- surv_first_flare_control_within$model.objects$models[["Exponential"]] # CS's choice is exponential

# Parametric survival model to extrapolate beyond trial (use both Weibull and exponential)
surv_first_flare_treat_beyond <- fit.fun(
  time = "time",
  status = "status",
  data = IPD_first_flare_data_treat,
  extrapolate = TRUE,
  times = seq(0, 1560, by = 52), 
  xlow = 0, xhigh = 1560,
  ylow = 0, yhigh = 1,
  mods = c("exp", "weibull")
)
weibull_model_treat_beyond <- surv_first_flare_treat_beyond$model.objects$models[["Weibull (AFT)"]]
exp_model_treat_beyond <- surv_first_flare_treat_beyond$model.objects$models[["Exponential"]]

surv_first_flare_control_beyond <- fit.fun(
  time = "time",
  status = "status",
  data = IPD_first_flare_data_control,
  extrapolate = TRUE,
  times = seq(0, 1560, by = 52),
  xlow = 52, xhigh = 1560,
  ylow = 0, yhigh = 1,
  mods = c("exp", "weibull")
)
weibull_model_control_beyond <- surv_first_flare_control_beyond$model.objects$models[["Weibull (AFT)"]]
exp_model_control_beyond <- surv_first_flare_control_beyond$model.objects$models[["Exponential"]]

# Function to calculate monthly transition probabilities from 0 to 30 years
calculate_transition_probs <- function(short_term_model, long_term_model) {
  
  # PART 1: Short-term monthly probabilities (0 to treatment duration)
  monthly_timepoints <- seq(0, 52, by = 52/12)
  
  # Get short-term survival probabilities
  short_surv_summary <- summary(short_term_model, type = "survival", t = monthly_timepoints, ci = TRUE)[[1]]
  short_surv_est <- short_surv_summary$est
  short_surv_lower <- short_surv_summary$lcl
  short_surv_upper <- short_surv_summary$ucl
  
  # Calculate monthly transition probabilities for first year
  n_short_months <- length(short_surv_est) - 1  # 12 months
  
  short_trans_prob_est <- 1 - (short_surv_est[-1] / short_surv_est[-length(short_surv_est)])
  short_trans_prob_lower <- 1 - (short_surv_upper[-1] / short_surv_upper[-length(short_surv_upper)])
  short_trans_prob_upper <- 1 - (short_surv_lower[-1] / short_surv_lower[-length(short_surv_lower)])
  
  # PART 2: Long-term yearly probabilities converted to monthly (years 2-30)
  yearly_timepoints <- seq(104, 1560, by = 52) # Years 2-30
  
  # Get long-term survival probabilities
  long_surv_summary <- summary(long_term_model, type = "survival", t = yearly_timepoints, ci = TRUE)[[1]]
  long_surv_est <- long_surv_summary$est
  long_surv_lower <- long_surv_summary$lcl
  long_surv_upper <- long_surv_summary$ucl
  
  # Add baseline from end of short-term period
  end_short_surv <- c(short_surv_est[length(short_surv_est)], long_surv_est)
  end_short_lower <- c(short_surv_lower[length(short_surv_lower)], long_surv_lower)
  end_short_upper <- c(short_surv_upper[length(short_surv_upper)], long_surv_upper)
  
  # Calculate yearly transition probabilities for long-term
  n_long_years <- length(end_short_surv) - 1  # 29 years (2-30)
  
  long_yearly_trans_est <- 1 - (end_short_surv[-1] / end_short_surv[-length(end_short_surv)])
  long_yearly_trans_lower <- 1 - (end_short_upper[-1] / end_short_upper[-length(end_short_upper)])
  long_yearly_trans_upper <- 1 - (end_short_lower[-1] / end_short_lower[-length(end_short_lower)])
  
  # Convert yearly to monthly equivalents
  long_monthly_trans_est <- 1 - (1 - long_yearly_trans_est)^(1/12)
  long_monthly_trans_lower <- 1 - (1 - long_yearly_trans_lower)^(1/12)
  long_monthly_trans_upper <- 1 - (1 - long_yearly_trans_upper)^(1/12)
  
  # PART 3: Combine into full 30-year monthly sequence
  
  # Create month sequence
  all_months <- 1:360  # 30 years * 12 months
  
  # Initialize results
  all_trans_prob <- numeric(360)
  all_trans_lower <- numeric(360)
  all_trans_upper <- numeric(360)
  
  # Fill in short-term (months 1-12)
  all_trans_prob[1:n_short_months] <- short_trans_prob_est
  all_trans_lower[1:n_short_months] <- short_trans_prob_lower
  all_trans_upper[1:n_short_months] <- short_trans_prob_upper
  
  # Fill in long-term (months 13-360)
  # Each yearly probability applies to 12 months
  for (i in 1:n_long_years) {
    start_month <- 12 + (i-1)*12 + 1
    end_month <- min(12 + i*12, 360)
    
    all_trans_prob[start_month:end_month] <- long_monthly_trans_est[i]
    all_trans_lower[start_month:end_month] <- long_monthly_trans_lower[i]
    all_trans_upper[start_month:end_month] <- long_monthly_trans_upper[i]
  }
  
  # Ensure valid probability bounds
  all_trans_lower <- pmax(0, all_trans_lower)
  all_trans_upper <- pmin(1, all_trans_upper)
  
  # Create final results dataframe
  results <- data.frame(
    month = all_months,
    year = ceiling(all_months / 12),
    trans_prob = all_trans_prob,
    trans_prob_lower = all_trans_lower,
    trans_prob_upper = all_trans_upper,
    period = c(rep("Short-term", 12), rep("Long-term", 348)),
    stringsAsFactors = FALSE
  )
  
  return(results)
}

# Transition probabolities calculation for each scenario
trans_probs_first_flare_treat_weibull <- calculate_transition_probs(
  exp_model_treat_within, 
  weibull_model_treat_beyond
)
trans_probs_first_flare_treat_exp <- calculate_transition_probs(
  exp_model_treat_within, 
  exp_model_treat_beyond
)
trans_probs_first_flare_control_weibull <- calculate_transition_probs(
  exp_model_control_within, 
  weibull_model_control_beyond
)
trans_probs_first_flare_control_exp <- calculate_transition_probs(
  exp_model_control_within, 
  exp_model_control_beyond
)

# Visualization function
plot_transition_probs <- function(trans_prob_df, title, max_months = 60) {
  df_subset <- trans_prob_df[1:min(max_months, nrow(trans_prob_df)), ]
  
  ggplot(df_subset, aes(x = month, y = trans_prob)) +
    geom_line(color = "blue", size = 1) +
    geom_ribbon(aes(ymin = trans_prob_lower, ymax = trans_prob_upper), 
                alpha = 0.3, fill = "blue") +
    labs(
      title = title,
      x = "Month",
      y = "Monthly Transition Probability"
    ) +
    theme_minimal() +
    scale_x_continuous(breaks = seq(0, max_months, 12)) +
    theme(plot.title = element_text(size = 12))
}

# Create plots
plot_transition_probs(trans_probs_first_flare_treat_weibull, 
                            "Treatment - Weibull extrapolation")
plot_transition_probs(trans_probs_first_flare_control_weibull, 
                            "Control - Weibull extrapolation")
plot_transition_probs(trans_probs_first_flare_treat_exp, 
                            "Treatment - Exponential extrapolation")
plot_transition_probs(trans_probs_first_flare_control_exp, 
                            "Control - Exponential extrapolation")


