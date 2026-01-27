##############################################
# Time to first flare transition probability #
##############################################

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

first_flare_newlydiag_data <- survival_digitize(img_path = here(getwd(),"01_data/Raw data/Time to first flare_newly diagnosed.png"), 
                                                num_curves = 4,
                                                x_start = 0,
                                                x_end = 52,
                                                x_increment = 4,
                                                y_start = 0,
                                                y_end = 100,
                                                y_increment = 20,
                                                y_text_vertical =  T)

first_flare_relapsed_data <- survival_digitize(img_path = here(getwd(),"01_data/Raw data/Time to first flare_relapsed.png"), 
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

first_flare_newlydiag_graph <- ggplot(data=first_flare_newlydiag_data) + geom_line(aes(x=time,y=St, colour=factor(curve))) + expand_limits(x = 0, y = 0) + scale_y_continuous(breaks=seq(0,100,20)) + scale_x_continuous(breaks=seq(0,52,4))
first_flare_newlydiag_graph

first_flare_relapsed_graph <- ggplot(data=first_flare_relapsed_data) + geom_line(aes(x=time,y=St, colour=factor(curve))) + expand_limits(x = 0, y = 0) + scale_y_continuous(breaks=seq(0,100,20)) + scale_x_continuous(breaks=seq(0,52,4))
first_flare_relapsed_graph

# Keeping the Km-curves generated for Tocilizumab QW + 26-week of tapering steroids and placebo + 52-week of tapering steroids
first_flare_data_itt_treat <- first_flare_itt_data[ which(first_flare_itt_data$curve==3),c(2,3) ] # treatment group
first_flare_data_itt_control <- first_flare_itt_data[ which(first_flare_itt_data$curve==1),c(2,3) ] # control group

first_flare_data_newlydiag_treat <- first_flare_newlydiag_data[ which(first_flare_newlydiag_data$curve==1),c(2,3) ] # treatment group
first_flare_data_newlydiag_control <- first_flare_newlydiag_data[ which(first_flare_newlydiag_data$curve==4),c(2,3) ] # control group

first_flare_data_relapsed_treat <- first_flare_relapsed_data[ which(first_flare_relapsed_data$curve==3),c(2,3) ] # treatment group
first_flare_data_relapsed_control <- first_flare_relapsed_data[ which(first_flare_relapsed_data$curve==4),c(2,3) ] # control group

# Generating the vector variables which provide information about the follow-up time and the number of patients at risk at every time point. 
weeks <- c(0,4,8,12,16,20,24,28,32,36,40,44,48,52)
n_risk_itt_treat <- c(100,93,88,85,85,81,77,74,71,69,67,64,63,5)
n_risk_itt_control <- c(51,48,44,41,38,35,32,30,28,25,22,17,15,0)
n_risk_newlydiag_treat <- c(47,43,42,40,40,40,38,36,35,34,33,33,32,4)
n_risk_newlydiag_control <- c(23,21,21,19,18,17,16,14,14,13,13,11,10,0)
n_risk_relapsed_treat <- c(53,50,46,45,45,41,39,38,36,35,34,31,31,1)
n_risk_relapsed_control <- c(28,27,23,22,20,18,16,16,14,12,9,6,5)

# Preprocess state of the KM curve information 
pre_first_flare_data_itt_treat <- preprocess(first_flare_data_itt_treat, trisk =weeks, nrisk=n_risk_itt_treat,totalpts = NULL,maxy=100)
pre_first_flare_data_itt_control <- preprocess(first_flare_data_itt_control, trisk =weeks, nrisk=n_risk_itt_control,totalpts = NULL,maxy=100)

pre_first_flare_data_newlydiag_treat <- preprocess(first_flare_data_newlydiag_treat, trisk =weeks, nrisk=n_risk_newlydiag_treat,totalpts = NULL,maxy=100)
pre_first_flare_data_newlydiag_control <- preprocess(first_flare_data_newlydiag_control, trisk =weeks, nrisk=n_risk_newlydiag_control,totalpts = NULL,maxy=100)

pre_first_flare_data_relapsed_treat <- preprocess(first_flare_data_relapsed_treat, trisk =weeks, nrisk=n_risk_relapsed_treat,totalpts = NULL,maxy=100)
pre_first_flare_data_relapsed_control <- preprocess(first_flare_data_relapsed_control, trisk =weeks, nrisk=n_risk_relapsed_control,totalpts = NULL,maxy=100)

# IPD generation. This dataset will be used for the time to event analysis with the flexsurv package
IPD_first_flare_data_itt_treat <- getIPD(prep=pre_first_flare_data_itt_treat,armID=1)$IPD
IPD_first_flare_data_itt_control <- getIPD(prep=pre_first_flare_data_itt_control,armID=0)$IPD

IPD_first_flare_data_newlydiag_treat <- getIPD(prep=pre_first_flare_data_newlydiag_treat,armID=1)$IPD
IPD_first_flare_data_newlydiag_control <- getIPD(prep=pre_first_flare_data_newlydiag_control,armID=0)$IPD

IPD_first_flare_data_relapsed_treat <- getIPD(prep=pre_first_flare_data_relapsed_treat,armID=1)$IPD
IPD_first_flare_data_relapsed_control <- getIPD(prep=pre_first_flare_data_relapsed_control,armID=0)$IPD

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

survfit2(Surv(time = time,  event = status)  ~ 1, data = IPD_first_flare_data_newlydiag_treat) %>%
  ggsurvfit() +
  labs(
    x = "Weeks",
    y = "Patients without flare -- Tocilizumab (ITT)"
  ) + 
  add_confidence_interval() +
  add_risktable()

summary(survfit(Surv(time, status) ~ 1, data = IPD_first_flare_data_newlydiag_treat), times = 52)

survfit2(Surv(time = time,  event = status)  ~ 1, data = IPD_first_flare_data_newlydiag_control) %>%
  ggsurvfit() +
  labs(
    x = "Weeks",
    y = "Patients without flare -- Steroid (ITT)"
  ) + 
  add_confidence_interval() +
  add_risktable()

summary(survfit(Surv(time, status) ~ 1, data = IPD_first_flare_data_newlydiag_control), times = 48.8)

survfit2(Surv(time = time,  event = status)  ~ 1, data = IPD_first_flare_data_relapsed_treat) %>%
  ggsurvfit() +
  labs(
    x = "Weeks",
    y = "Patients without flare -- Tocilizumab (Relapsed)"
  ) + 
  add_confidence_interval() +
  add_risktable()

summary(survfit(Surv(time, status) ~ 1, data = IPD_first_flare_data_relapsed_treat), times = 52)

survfit2(Surv(time = time,  event = status)  ~ 1, data = IPD_first_flare_data_relapsed_control) %>%
  ggsurvfit() +
  labs(
    x = "Weeks",
    y = "Patients without flare -- Steroid (Relapsed)"
  ) + 
  add_confidence_interval() +
  add_risktable()

summary(survfit(Surv(time, status) ~ 1, data = IPD_first_flare_data_relapsed_control), times = 52)


# Fit multiple survival models on within-trial data
surv_first_flare_itt_treat_within <- fit.fun(
  time = "time",
  status = "status",
  data = IPD_first_flare_data_itt_treat,
  xlow = 0, xhigh = 52,
  ylow = 0.75, yhigh = 1,
  mods = c("exp", "weibull", "gamma", "lnorm")
)
surv_first_flare_itt_treat_within$AIC # Although exponential has the smallest AIC, Weibull fits better visually.
weibull_model_itt_treat_within <- surv_first_flare_itt_treat_within$model.objects$models[["Weibull (AFT)"]] # CS's choice

surv_first_flare_itt_control_within <- fit.fun(
  time = "time", 
  status = "status", 
  data = IPD_first_flare_data_itt_control,
  xlow = 0, xhigh = 52,
  ylow = 0.30, yhigh = 1,
  mods = c("exp", "weibull", "gamma", "lnorm")
) 
surv_first_flare_itt_control_within$AIC # Log-normal (but exponential is the 2nd smallest)
exp_model_itt_control_within <- surv_first_flare_itt_control_within$model.objects$models[["Exponential"]] # CS's choice 

surv_first_flare_newlydiag_treat_within <- fit.fun(
  time = "time",
  status = "status",
  data = IPD_first_flare_data_newlydiag_treat,
  xlow = 0, xhigh = 52,
  ylow = 0.70, yhigh = 1,
  mods = c("exp", "weibull", "gamma", "lnorm")
)
surv_first_flare_newlydiag_treat_within$AIC # Exponential is the largest, but similar in magnitudes
exp_model_newlydiag_treat_within <- surv_first_flare_newlydiag_treat_within$model.objects$models[["Exponential"]] # CS's choice 

surv_first_flare_newlydiag_control_within <- fit.fun(
  time = "time", 
  status = "status", 
  data = IPD_first_flare_data_newlydiag_control,
  xlow = 0, xhigh = 52,
  ylow = 0.40, yhigh = 1,
  mods = c("exp", "weibull", "gamma", "lnorm")
) 
surv_first_flare_newlydiag_control_within$AIC # Exponential
exp_model_newlydiag_control_within <- surv_first_flare_newlydiag_control_within$model.objects$models[["Exponential"]] # CS's choice 

surv_first_flare_relapsed_treat_within <- fit.fun(
  time = "time",
  status = "status",
  data = IPD_first_flare_data_relapsed_treat,
  xlow = 0, xhigh = 52,
  ylow = 0.75, yhigh = 1,
  mods = c("exp", "weibull", "gamma", "lnorm")
)
surv_first_flare_relapsed_treat_within$AIC # Exponential 
exp_model_relapsed_treat_within <- surv_first_flare_relapsed_treat_within$model.objects$models[["Exponential"]] # CS's choice

surv_first_flare_relapsed_control_within <- fit.fun(
  time = "time", 
  status = "status", 
  data = IPD_first_flare_data_relapsed_control,
  xlow = 0, xhigh = 52,
  ylow = 0.30, yhigh = 1,
  mods = c("exp", "weibull", "gamma", "lnorm")
) 
surv_first_flare_relapsed_control_within$AIC # Although Weibull and Gamma are smaller than log-normal, CS selected log-normal.
lnorm_model_relapsed_control_within <- surv_first_flare_relapsed_control_within$model.objects$models[["log-Normal"]] # CS's choice 

# Parametric survival model to extrapolate beyond trial (weibull)
surv_first_flare_itt_treat_beyond <- fit.fun(
  time = "time",
  status = "status",
  data = IPD_first_flare_data_itt_treat,
  extrapolate = TRUE,
  times = seq(0, time_horizon_weeks, by = 52), 
  xlow = 0, xhigh = time_horizon_weeks,
  ylow = 0, yhigh = 1,
  mods = c("exp", "weibull")
)
weibull_model_itt_treat_beyond <- surv_first_flare_itt_treat_beyond$model.objects$models[["Weibull (AFT)"]]

surv_first_flare_itt_control_beyond <- fit.fun(
  time = "time",
  status = "status",
  data = IPD_first_flare_data_itt_control,
  extrapolate = TRUE,
  times = seq(0, time_horizon_weeks, by = 52),
  xlow = 0, xhigh = time_horizon_weeks,
  ylow = 0, yhigh = 1,
  mods = c("exp", "weibull")
)
weibull_model_itt_control_beyond <- surv_first_flare_itt_control_beyond$model.objects$models[["Weibull (AFT)"]]

surv_first_flare_newlydiag_treat_beyond <- fit.fun(
  time = "time",
  status = "status",
  data = IPD_first_flare_data_newlydiag_treat,
  extrapolate = TRUE,
  times = seq(0, time_horizon_weeks, by = 52), 
  xlow = 0, xhigh = time_horizon_weeks,
  ylow = 0, yhigh = 1,
  mods = c("exp", "weibull")
)
weibull_model_newlydiag_treat_beyond <- surv_first_flare_newlydiag_treat_beyond$model.objects$models[["Weibull (AFT)"]]

surv_first_flare_newlydiag_control_beyond <- fit.fun(
  time = "time",
  status = "status",
  data = IPD_first_flare_data_newlydiag_control,
  extrapolate = TRUE,
  times = seq(0, time_horizon_weeks, by = 52),
  xlow = 0, xhigh = time_horizon_weeks,
  ylow = 0, yhigh = 1,
  mods = c("exp", "weibull")
)
weibull_model_newlydiag_control_beyond <- surv_first_flare_newlydiag_control_beyond$model.objects$models[["Weibull (AFT)"]]

surv_first_flare_relapsed_treat_beyond <- fit.fun(
  time = "time",
  status = "status",
  data = IPD_first_flare_data_relapsed_treat,
  extrapolate = TRUE,
  times = seq(0, time_horizon_weeks, by = 52), 
  xlow = 0, xhigh = time_horizon_weeks,
  ylow = 0, yhigh = 1,
  mods = c("exp", "weibull", "lnorm")
)
weibull_model_relapsed_treat_beyond <- surv_first_flare_relapsed_treat_beyond$model.objects$models[["Weibull (AFT)"]]

surv_first_flare_relapsed_control_beyond <- fit.fun(
  time = "time",
  status = "status",
  data = IPD_first_flare_data_relapsed_control,
  extrapolate = TRUE,
  times = seq(0, time_horizon_weeks, by = 52),
  xlow = 0, xhigh = time_horizon_weeks,
  ylow = 0, yhigh = 1,
  mods = c("exp", "weibull", "lnorm")
)
weibull_model_relapsed_control_beyond <- surv_first_flare_relapsed_control_beyond$model.objects$models[["Weibull (AFT)"]]

# Function to calculate monthly transition probabilities 
calculate_transition_probs <- function(short_term_model, long_term_model) {
  
  # Short-term monthly probabilities (0 to treatment duration)
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
  
  # 2. Long-term yearly probabilities converted to monthly (beyond treatment duration)
  yearly_timepoints <- seq(104, time_horizon_weeks, by = 52) # Years 2-10
  
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
  n_long_years <- length(end_short_surv) - 1  
  
  long_yearly_trans_est <- 1 - (end_short_surv[-1] / end_short_surv[-length(end_short_surv)])
  long_yearly_trans_lower <- 1 - (end_short_upper[-1] / end_short_upper[-length(end_short_upper)])
  long_yearly_trans_upper <- 1 - (end_short_lower[-1] / end_short_lower[-length(end_short_lower)])
  
  # Convert yearly to monthly equivalents
  long_monthly_trans_est <- 1 - (1 - long_yearly_trans_est)^(1/12)
  long_monthly_trans_lower <- 1 - (1 - long_yearly_trans_lower)^(1/12)
  long_monthly_trans_upper <- 1 - (1 - long_yearly_trans_upper)^(1/12)
  
  # 3. Combine the short-term and long-term monthly sequence
  
  # Create month sequence
  all_months <- 1:n_cycles  
  
  # Initialize results
  all_trans_prob <- numeric(n_cycles)
  all_trans_lower <- numeric(n_cycles)
  all_trans_upper <- numeric(n_cycles)
  
  # Fill in short-term (during the treatment duration)
  all_trans_prob[1:n_short_months] <- short_trans_prob_est
  all_trans_lower[1:n_short_months] <- short_trans_prob_lower
  all_trans_upper[1:n_short_months] <- short_trans_prob_upper
  
  # Fill in long-term (beyond the treatment duration)
  # Each yearly probability applies to 12 months
  for (i in 1:n_long_years) {
    start_month <- 12 + (i-1)*12 + 1
    end_month <- min(12 + i*12, n_cycles)
    
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
    period = c(rep("Short-term", 12), rep("Long-term", n_cycles-12)),
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

trans_probs_first_flare_newlydiag_treat <- calculate_transition_probs(
  exp_model_newlydiag_treat_within, 
  weibull_model_newlydiag_treat_beyond
)
trans_probs_first_flare_newlydiag_control <- calculate_transition_probs(
  exp_model_newlydiag_control_within, 
  weibull_model_newlydiag_control_beyond
)

trans_probs_first_flare_relapsed_treat <- calculate_transition_probs(
  exp_model_relapsed_treat_within, 
  weibull_model_relapsed_treat_beyond
)
trans_probs_first_flare_relapsed_control <- calculate_transition_probs(
  lnorm_model_relapsed_control_within, 
  weibull_model_relapsed_control_beyond
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
plot_transition_probs(trans_probs_first_flare_newlydiag_treat, 
                      "Tocilizumab Newly Diagnosed")
plot_transition_probs(trans_probs_first_flare_newlydiag_control, 
                      "Placebo Newly Diagnosed")
plot_transition_probs(trans_probs_first_flare_relapsed_treat, 
                      "Tocilizumab Relapsed")
plot_transition_probs(trans_probs_first_flare_relapsed_control, 
                      "Placebo Relapsed")


