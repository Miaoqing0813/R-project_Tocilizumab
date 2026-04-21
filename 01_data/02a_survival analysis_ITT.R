##############################################################
# Time to first flare transition probability (ITT population)#
##############################################################

time_horizon_weeks <- 5*52
time_horizon_months <- 5*12
treatment_duration_weeks <- 52

# To digitize from published Kaplan–Meier (KM) curves
first_flare_itt_data <- survival_digitize(img_path = here(getwd(),"01_data/Raw data/Time to first flare_ITT.png"), 
                                      num_curves = 4,
                                      x_start = 0,
                                      x_end = 52,
                                      x_increment = 4,
                                      y_start = 0,
                                      y_end = 100,
                                      y_increment = 20,
                                      y_text_vertical =  T)
                                      
# Generating a graph and checking how it looks compared to the original
first_flare_itt_graph <- ggplot(data=first_flare_itt_data) + geom_line(aes(x=time,y=St, colour=factor(curve))) + expand_limits(x = 0, y = 0) + scale_y_continuous(breaks=seq(0,100,5)) + scale_x_continuous(breaks=seq(0,132,6))
first_flare_itt_graph

# Keeping the Km-curves generated for Tocilizumab QW + 26-week of tapering steroids and placebo + 52-week of tapering steroids
first_flare_data_itt_treat <- first_flare_itt_data[ which(first_flare_itt_data$curve==3),c(2,3) ] # treatment group
first_flare_data_itt_control <- first_flare_itt_data[ which(first_flare_itt_data$curve==1),c(2,3) ] # control group

# Generating the vector variables which provide information about the follow-up time and the number of patients at risk at every time point. 
weeks <- c(0,4,8,12,16,20,24,28,32,36,40,44,48,52)
n_risk_itt_treat <- c(100,93,88,85,85,81,77,74,71,69,67,64,63,5)
n_risk_itt_control <- c(51,48,44,41,38,35,32,30,28,25,22,17,15,0)

#Preprocess state of the KM curve information 
pre_first_flare_data_itt_treat <- preprocess(first_flare_data_itt_treat, trisk =weeks, nrisk=n_risk_itt_treat,totalpts = NULL,maxy=100)
pre_first_flare_data_itt_control <- preprocess(first_flare_data_itt_control, trisk =weeks, nrisk=n_risk_itt_control,totalpts = NULL,maxy=100)

# IPD generation. This dataset will be used for the time to event analysis with the flexsurv package
IPD_first_flare_data_itt_treat <- getIPD(prep=pre_first_flare_data_itt_treat,armID=1)$IPD
IPD_first_flare_data_itt_control <- getIPD(prep=pre_first_flare_data_itt_control,armID=0)$IPD
IPD_combined <- rbind(
  data.frame(IPD_first_flare_data_itt_treat,   arm = 1),
  data.frame(IPD_first_flare_data_itt_control, arm = 0)
)
IPD_combined$arm <- factor(IPD_combined$arm, levels = c(0, 1))

# KM-curve -- within trial period
survfit2(Surv(time = time,  event = status)  ~ 1, data = IPD_first_flare_data_itt_treat) %>%
  ggsurvfit() +
  labs(
    x = "Weeks",
    y = "Patients without flare -- Tocilizumab (ITT)"
  ) + 
  add_confidence_interval() +
  add_risktable()

summary(survfit(Surv(time, status) ~ 1, data = IPD_first_flare_data_itt_treat), times = 52)

survfit2(Surv(time = time,  event = status)  ~ 1, data = IPD_first_flare_data_itt_control) %>%
  ggsurvfit() +
  labs(
    x = "Weeks",
    y = "Patients without flare -- Steroid (ITT)"
  ) + 
  add_confidence_interval() +
  add_risktable()

summary(survfit(Surv(time, status) ~ 1, data = IPD_first_flare_data_itt_control), times = 51.7)

survfit2(Surv(time = time,  event = status)  ~ arm, data = IPD_combined) %>%
  ggsurvfit() +
  labs(
    x = "Weeks",
    y = "Patients without flare -- Steroid (ITT)"
  ) + 
  add_confidence_interval() 

# Fit multiple survival models during treatment period
surv_first_flare_itt_treat_within <- fit.fun(
  time = "time",
  status = "status",
  data = IPD_first_flare_data_itt_treat,
  xlow = 0, xhigh = treatment_duration_weeks,
  ylow = 0.75, yhigh = 1,
  mods = c("exp", "weibull", "gamma", "lnorm")
)
surv_first_flare_itt_treat_within$AIC # Although exponential has the smallest AIC, Weibull fits better visually.
weibull_model_itt_treat_within <- surv_first_flare_itt_treat_within$model.objects$models[["Weibull (AFT)"]] # CS's choice

surv_first_flare_itt_control_within <- fit.fun(
  time = "time", 
  status = "status", 
  data = IPD_first_flare_data_itt_control,
  xlow = 0, xhigh = treatment_duration_weeks,
  ylow = 0.30, yhigh = 1,
  mods = c("exp", "weibull", "gamma", "lnorm")
) 
surv_first_flare_itt_control_within$AIC # Log-normal (but exponential is the 2nd smallest)
exp_model_itt_control_within <- surv_first_flare_itt_control_within$model.objects$models[["Exponential"]] # CS's choice 

# Parametric survival model to extrapolate beyond trial (weibull)
surv_first_flare_itt_treat_beyond <- fit.fun(
  time = "time",
  status = "status",
  data = IPD_first_flare_data_itt_treat,
  extrapolate = TRUE,
  times = seq(0, time_horizon_weeks, by = 52), 
  xlow = 0, xhigh = time_horizon_weeks,
  ylow = 0, yhigh = 1,
  mods = c("exp", "weibull", "gamma", "lnorm")
)
weibull_model_itt_treat_beyond <- surv_first_flare_itt_treat_beyond$model.objects$models[["Weibull (AFT)"]]
exp_model_itt_treat_beyond <- surv_first_flare_itt_treat_beyond$model.objects$models[["Exponential"]]
summary(weibull_model_itt_treat_beyond, t = 152, type = "survival") # compare with the 3-year GiACTA trial results

surv_first_flare_itt_control_beyond <- fit.fun(
  time = "time",
  status = "status",
  data = IPD_first_flare_data_itt_control,
  extrapolate = TRUE,
  times = seq(0, time_horizon_weeks, by = 52),
  xlow = 0, xhigh = time_horizon_weeks,
  ylow = 0, yhigh = 1,
  mods = c("exp", "weibull", "gamma", "lnorm")
)
weibull_model_itt_control_beyond <- surv_first_flare_itt_control_beyond$model.objects$models[["Weibull (AFT)"]]
exp_model_itt_control_beyond <- surv_first_flare_itt_control_beyond$model.objects$models[["Exponential"]]
summary(weibull_model_itt_control_beyond, t = 152, type = "survival") # compare with the 3-year GiACTA trial results

# Function to calculate monthly transition probabilities
calculate_transition_probs <- function(short_term_model, long_term_model) {
  
  # Short-term monthly probabilities (0 to treatment duration)
  monthly_timepoints <- seq(0, treatment_duration_weeks, by = 52/12)
  
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
  
  # 2. Long-term yearly probabilities converted to monthly (beyond treatment duration)
  yearly_timepoints <- seq(treatment_duration_weeks+52, time_horizon_weeks, by = 52) 
  
  # Get long-term survival probabilities
  long_surv_summary <- summary(long_term_model, type = "survival", t = yearly_timepoints, ci = TRUE)[[1]]
  long_surv_est <- long_surv_summary$est
  long_surv_lower <- long_surv_summary$lcl
  long_surv_upper <- long_surv_summary$ucl
  
  # Add baseline from end of short-term period
  endshort_long_surv_est <- c(short_surv_est[length(short_surv_est)], long_surv_est)
  endshort_long_surv_lower <- c(short_surv_lower[length(short_surv_lower)], long_surv_lower)
  endshort_long_surv_upper <- c(short_surv_upper[length(short_surv_upper)], long_surv_upper)
  
  # Calculate yearly transition probabilities for long-term
  n_long_years <- length(endshort_long_surv_est) - 1  
  
  long_yearly_trans_est <- 1 - (endshort_long_surv_est[-1] / endshort_long_surv_est[-length(endshort_long_surv_est)])
  long_yearly_trans_lower <- 1 - (endshort_long_surv_lower[-1] / endshort_long_surv_lower[-length(endshort_long_surv_lower)])
  long_yearly_trans_upper <- 1 - (endshort_long_surv_upper[-1] / endshort_long_surv_upper[-length(endshort_long_surv_upper)])
  
  # Convert yearly to monthly equivalents
  long_monthly_trans_est <- 1 - (1 - long_yearly_trans_est)^(1/12)
  long_monthly_trans_lower <- 1 - (1 - long_yearly_trans_lower)^(1/12)
  long_monthly_trans_upper <- 1 - (1 - long_yearly_trans_upper)^(1/12)
  
  # 3. Combine the short-term and long-term monthly sequence
  
  # Create month sequence
  all_months <- 1:time_horizon_months  
  
  # Initialize results
  all_trans_prob <- numeric(time_horizon_months)
  all_trans_lower <- numeric(time_horizon_months)
  all_trans_upper <- numeric(time_horizon_months)
  
  # Fill in short-term (during the treatment duration)
  all_trans_prob[1:n_short_months] <- short_trans_prob_est
  all_trans_lower[1:n_short_months] <- short_trans_prob_lower
  all_trans_upper[1:n_short_months] <- short_trans_prob_upper
  
  # Fill in long-term (beyond the treatment duration)
  # Each yearly probability applies to 12 months
  for (i in 1:n_long_years) {
    start_month <- 12 + (i-1)*12 + 1
    end_month <- min(12 + i*12, time_horizon_months)
    
    all_trans_prob[start_month:end_month] <- long_monthly_trans_est[i]
    all_trans_lower[start_month:end_month] <- long_monthly_trans_lower[i]
    all_trans_upper[start_month:end_month] <- long_monthly_trans_upper[i]
  }
  
  # Ensure valid probability bounds
  all_trans_lower <- pmax(0, all_trans_lower)
  all_trans_upper <- pmin(1, all_trans_upper)
  
  # Create final results dataframe
  trans_prob <- data.frame(
    month = all_months,
    year = ceiling(all_months / 12),
    trans_prob = all_trans_prob,
    trans_prob_lower = all_trans_lower,
    trans_prob_upper = all_trans_upper,
    period = c(rep("Short-term", 12), rep("Long-term", time_horizon_months-12)),
    stringsAsFactors = FALSE
  )
  
  return(trans_prob)
}

# Transition probabilities calculation for each scenario
trans_probs_first_flare_itt_treat <- calculate_transition_probs(
  weibull_model_itt_treat_within, 
  weibull_model_itt_treat_beyond
)
trans_probs_first_flare_itt_control <- calculate_transition_probs(
  exp_model_itt_control_within, 
  weibull_model_itt_control_beyond
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
plot_transition_probs(trans_probs_first_flare_itt_treat, 
                      "Tocilizumab ITT")
plot_transition_probs(trans_probs_first_flare_itt_control, 
                      "Placebo ITT")


