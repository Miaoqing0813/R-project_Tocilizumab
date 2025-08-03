###############################################
# Time to first flare transition probability #
###############################################

# To digitize from published Kaplan–Meier (KM) curves
first_flare_data <- survival_digitize(img_path = here(getwd(),"01_data/Raw data/Remission to first flare_all patients_clean.png"), 
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
first_flare_graph <- ggplot(data=first_flare_data) + geom_line(aes(x=time,y=St, colour=factor(curve))) + expand_limits(x = 0, y = 0) + scale_y_continuous(breaks=seq(0,100,5)) + scale_x_continuous(breaks=seq(0,132,6))
first_flare_graph

# Keeping the Km-curves generated for Tocilizumab QW + 26-week of tapering steroids and placebo + 52-week of tapering steroids
first_flare_data_treat <- first_flare_data[ which(first_flare_data$curve==3),c(2,3) ] # treatment group
first_flare_data_control <- first_flare_data[ which(first_flare_data$curve==1),c(2,3) ] # control group

# Generating the vector variables which provide information about the follow-up time 
# and the number of patients at risk at every time point. Both of these vectors
# will be used for the preprocess state for generating IPD datasets.
weeks <- c(0,4,8,12,16,20,24,28,32,36,40,44,48,52)
n_risk_treat <- c(100,93,88,85,85,81,77,74,71,69,67,64,63,5)
n_risk_control <- c(51,48,44,41,38,35,32,30,28,25,22,17,15,0)

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
summary(survfit(Surv(time, status) ~ 1, data = IPD_first_flare_data_control), times = 51.5)

# Fit multiple survival models on within-trial data
surv_first_flare_treat_within <- fit.fun(
  time = "time",
  status = "status",
  data = IPD_first_flare_data_treat,
  xlow = 0, xhigh = 52,
  ylow = 0.75, yhigh = 1,
  mods = c("exp", "weibull", "gamma", "lnorm")
)
surv_first_flare_treat_within$AIC # Exponential
surv_first_flare_treat_within$BIC # Exponential
weibull_model_treat_within <- surv_first_flare_treat_within$model.objects$models[["Weibull (AFT)"]] # from the graph, Weibull actually fits better; and this is the CS's choice as well.

surv_first_flare_control_within <- fit.fun(
  time = "time", 
  status = "status", 
  data = IPD_first_flare_data_control,
  xlow = 0, xhigh = 52,
  ylow = 0.30, yhigh = 1,
  mods = c("exp", "weibull", "gamma", "lnorm")
) 
surv_first_flare_control_within$AIC # Log-normal (but exponential is the 2nd smallest)
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

# Let's apply different survival models for the treatment period and the extrapolation period
## Assume 1-year treatment duration
monthly_within_1yr <- seq(0, 52, by = 52/12) # monthly from week 0 to 52
# Use observed KM-derived survival at 0 and 52 weeks 
km_treat_summary <- summary(km_treat, times = monthly_within_1yr, extend = TRUE)
first_flare_treat_within_1yr <- data.frame(
  time = seq(0, 12),
  surv = km_treat_summary$surv,
  lower = km_treat_summary$lower,
  upper = km_treat_summary$upper
)
km_control_summary <- summary(km_control, times = monthly_within_1yr, extend = TRUE)
first_flare_control_within_1yr <- data.frame(
  time = seq(0, 12),
  surv = km_control_summary$surv,
  lower = km_control_summary$lower,
  upper = km_control_summary$upper
)
# Extrapolate using Weibull and exponential models
monthly_beyond_1yr <- seq(52+52/12, 1560, by = 52/12)  # monthly from year 2 to year 30
weibull_treat_summary <- summary(weibull_model_treat_beyond, type = "survival", t = monthly_beyond_1yr, ci = TRUE)[[1]]
first_flare_treat_beyond_1yr_weibull <- data.frame(
  time = seq(13, 360),
  surv = weibull_treat_summary$est,
  lower = weibull_treat_summary$lcl,
  upper = weibull_treat_summary$ucl
)
weibull_control_summary <- summary(weibull_model_control_beyond, type = "survival", t = monthly_beyond_1yr, ci = TRUE)[[1]]
first_flare_control_beyond_1yr_weibull <- data.frame(
  time = seq(13, 360),
  surv = weibull_control_summary$est,
  lower = weibull_control_summary$lcl,
  upper = weibull_control_summary$ucl
)
exp_treat_summary <- summary(exp_model_treat_beyond, type = "survival", t = monthly_beyond_1yr, ci = TRUE)[[1]]
first_flare_treat_beyond_1yr_exp <- data.frame(
  time = seq(13, 360),
  surv = exp_treat_summary$est,
  lower = exp_treat_summary$lcl,
  upper = exp_treat_summary$ucl
)
exp_control_summary <- summary(exp_model_control_beyond, type = "survival", t = monthly_beyond_1yr, ci = TRUE)[[1]]
first_flare_control_beyond_1yr_exp <- data.frame(
  time = seq(13, 360),
  surv = exp_control_summary$est,
  lower = exp_control_summary$lcl,
  upper = exp_control_summary$ucl
)
# Combine
surv_probs_first_flare_treat_1yr_weibull <- rbind(first_flare_treat_within_1yr, first_flare_treat_beyond_1yr_weibull)
surv_probs_first_flare_treat_1yr_exp <- rbind(first_flare_treat_within_1yr, first_flare_treat_beyond_1yr_exp)
surv_probs_first_flare_control_1yr_weibull <- rbind(first_flare_control_within_1yr, first_flare_control_beyond_1yr_weibull)
surv_probs_first_flare_control_1yr_exp <- rbind(first_flare_control_within_1yr, first_flare_control_beyond_1yr_exp)

## Assume 2-year treatment duration
monthly_2yr <- seq(0, 104, by = 52/12)
# Weibull-based survival for the treatment group and exponential for the control group
weibull_treat_summary <- summary(weibull_model_treat_within, type = "survival", t = monthly_2yr, ci = TRUE)[[1]]
first_flare_treat_within_2yr <- data.frame(
  time = seq(0, 24),
  surv = weibull_treat_summary$est,
  lower = weibull_treat_summary$lcl,
  upper = weibull_treat_summary$ucl
)
exp_control_summary <- summary(exp_model_control_within, type = "survival", t = monthly_2yr, ci = TRUE)[[1]]
first_flare_control_within_2yr <- data.frame(
  time = seq(0, 24),
  surv = weibull_treat_summary$est,
  lower = weibull_treat_summary$lcl,
  upper = weibull_treat_summary$ucl
)
# Extrapolate using Weibull and exponential models
monthly_beyond_2yr <- seq(104+52/12, 1560, by = 52/12)  # monthly from year 3 to year 30

weibull_treat_summary <- summary(weibull_model_treat_beyond, type = "survival", t = monthly_beyond_2yr, ci = TRUE)[[1]]
first_flare_treat_beyond_2yr_weibull <- data.frame(
  time = seq(25, 360),
  surv = weibull_treat_summary$est,
  lower = weibull_treat_summary$lcl,
  upper = weibull_treat_summary$ucl
)
weibull_control_summary <- summary(weibull_model_control_beyond, type = "survival", t = monthly_beyond_2yr, ci = TRUE)[[1]]
first_flare_control_beyond_2yr_weibull <- data.frame(
  time = seq(25, 360),
  surv = weibull_control_summary$est,
  lower = weibull_control_summary$lcl,
  upper = weibull_control_summary$ucl
)
exp_treat_summary <- summary(exp_model_treat_beyond, type = "survival", t = monthly_beyond_2yr, ci = TRUE)[[1]]
first_flare_treat_beyond_2yr_exp <- data.frame(
  time = seq(25, 360),
  surv = exp_treat_summary$est,
  lower = exp_treat_summary$lcl,
  upper = exp_treat_summary$ucl
)
exp_control_summary <- summary(exp_model_control_beyond, type = "survival", t = monthly_beyond_2yr, ci = TRUE)[[1]]
first_flare_control_beyond_2yr_exp <- data.frame(
  time = seq(25, 360),
  surv = exp_control_summary$est,
  lower = exp_control_summary$lcl,
  upper = exp_control_summary$ucl
)
# Combine
surv_probs_first_flare_treat_2yr_weibull <- rbind(first_flare_treat_within_2yr, first_flare_treat_beyond_2yr_weibull)
surv_probs_first_flare_control_2yr_weibull <- rbind(first_flare_control_within_2yr, first_flare_control_beyond_2yr_weibull)
surv_probs_first_flare_treat_2yr_exp <- rbind(first_flare_treat_within_2yr, first_flare_treat_beyond_2yr_exp)
surv_probs_first_flare_control_2yr_exp <- rbind(first_flare_control_within_2yr, first_flare_control_beyond_2yr_exp)

# Calculate transition probabilities
trans_prob_df <- function(surv_df, surv_floor = 0.0001) {
  
  surv <- surv_df$surv
  lower <- surv_df$lower
  upper <- surv_df$upper
  
  # Truncate survival and CI values below floor to zero
  surv[surv < surv_floor] <- 0
  lower[lower < surv_floor] <- 0
  upper[upper < surv_floor] <- 0
  
  # Compute transition probabilities
  tp <- 1 - (surv[-1] / surv[-length(surv)])
  tp_lower <- 1 - (upper[-1] / upper[-length(upper)])
  tp_upper <- 1 - (lower[-1] / lower[-length(lower)])
  
  # Fix function for NA and logical cases
  fix_tp <- function(x, S_t, S_tp1) {
    x[is.na(x) & S_t == 0 & S_tp1 == 0] <- 1
    
    return(x)
  }
  
  # Clean transition estimates
  tp <- fix_tp(tp, surv[-length(surv)], surv[-1])
  tp_lower <- fix_tp(tp_lower, upper[-length(upper)], upper[-1])
  tp_upper <- fix_tp(tp_upper, lower[-length(lower)], lower[-1])
  
  # If tp == 1, set bounds to 1 as well
  tp_is_1 <- tp == 1
  tp_lower[tp_is_1] <- 1
  tp_upper[tp_is_1] <- 1
  
  # Return data frame
  data.frame(
    interval = 1:(nrow(surv_df) - 1),
    trans_prob = tp,
    trans_prob_lower = tp_lower,
    trans_prob_upper = tp_upper
  )
}

trans_probs_first_flare_treat_1yr_weibull <- trans_prob_df(surv_probs_first_flare_treat_1yr_weibull)
trans_probs_first_flare_control_1yr_weibull <- trans_prob_df(surv_probs_first_flare_control_1yr_weibull)
trans_probs_first_flare_treat_1yr_exp <- trans_prob_df(surv_probs_first_flare_treat_1yr_exp)
trans_probs_first_flare_control_1yr_exp <- trans_prob_df(surv_probs_first_flare_control_1yr_exp)

trans_probs_first_flare_treat_2yr_weibull <- trans_prob_df(surv_probs_first_flare_treat_2yr_weibull)
trans_probs_first_flare_control_2yr_weibull <- trans_prob_df(surv_probs_first_flare_control_2yr_weibull)
trans_probs_first_flare_treat_2yr_exp <- trans_prob_df(surv_probs_first_flare_treat_2yr_exp)
trans_probs_first_flare_control_2yr_exp <- trans_prob_df(surv_probs_first_flare_control_2yr_exp)




