###############################################
# Time to first flare transition probability #
###############################################

# codes taken from https://rdrr.io/github/DARTH-git/darthtools/src/R/survival_functions.R

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

#Generating a graph and checking how it looks compared to the original
first_flare_graph <- ggplot(data=first_flare_data) + geom_line(aes(x=time,y=St, colour=factor(curve))) + expand_limits(x = 0, y = 0) + scale_y_continuous(breaks=seq(0,100,5)) + scale_x_continuous(breaks=seq(0,132,6))
first_flare_graph

# Keeping the Km-curves generated for Tocilizumab QW + 26-week of tapering steroids and placebo + 52-week of tapering steroids
first_flare_data_treat <- first_flare_data[ which(first_flare_data$curve==3),c(2,3) ] # treatment group
first_flare_data_control <- first_flare_data[ which(first_flare_data$curve==1),c(2,3) ] # control group

#Generating the vector variables which provide information about the follow-up time 
#and the number of patients at risk at every time point. Both of these vectors
#will be used for the preprocess state for generating IPD datasets.
first_flare_weeks <- c(0,4,8,12,16,20,24,28,32,36,40,44,48,52)
first_flare_n_risk_treat <- c(100,93,88,85,85,81,77,74,71,69,67,64,63,5)
first_flare_n_risk_control <- c(51,48,44,41,38,35,32,30,28,25,22,17,15,0)

#Preprocess state of the KM curve information 
pre_first_flare_data_treat <- preprocess(first_flare_data_treat, trisk =first_flare_weeks, nrisk=first_flare_n_risk_treat,totalpts = NULL,maxy=100)
pre_first_flare_data_control <- preprocess(first_flare_data_control, trisk =first_flare_weeks, nrisk=first_flare_n_risk_control,totalpts = NULL,maxy=100)

# Dataset IPD generation. This dataset will be used for the time to event analysis with the flexsurv package
est_first_flare_data_treat <- getIPD(prep=pre_first_flare_data_treat,armID=1,tot.events=NULL)
est_first_flare_data_control <- getIPD(prep=pre_first_flare_data_control,armID=0,tot.events=NULL)
IPD_first_flare_data_treat <- est_first_flare_data_treat$IPD
IPD_first_flare_data_control <- est_first_flare_data_control$IPD
names(IPD_first_flare_data_treat) <- c("time","status","treat")
names(IPD_first_flare_data_control) <- c("time","status","treat")

# KM-curve
survfit2(Surv(time = time,  event = status)  ~ 1, data = IPD_first_flare_data_treat) %>%
  ggsurvfit() +
  labs(
    x = "Weeks",
    y = "Overall survival probability for the Tocilizumab group"
  ) + 
  add_confidence_interval() +
  add_risktable()

survfit2(Surv(time = time,  event = status)  ~ 1, data = IPD_first_flare_data_control) %>%
  ggsurvfit() +
  labs(
    x = "Weeks",
    y = "Overall survival probability for the comparison group"
  ) + 
  add_confidence_interval() +
  add_risktable()

# Estimating x-year survival and median survival time
km_treat <- survfit(Surv(time, status) ~ 1, data = IPD_first_flare_data_treat) 
km_control <- survfit(Surv(time, status) ~ 1, data = IPD_first_flare_data_control)  
summary(survfit(Surv(time, status) ~ 1, data = IPD_first_flare_data_treat), times = 52)
summary(survfit(Surv(time, status) ~ 1, data = IPD_first_flare_data_control), times = 51)

# Fit multiple survival models on survival data
surv_first_flare_treat_original <- fit.fun(time = "time", status = "status", data = IPD_first_flare_data_treat,
                                     xlow = 0, xhigh = 53,
                                     ylow = 0.75, yhigh = 1,
                                     mods = c("exp", "weibull", "gamma", "lnorm")
                                     ) 

surv_first_flare_treat <- fit.fun( time = "time", status = "status", data = IPD_first_flare_data_treat,
                                   extrapolate = TRUE, times=seq(0, 1560, by = 100),
                                   xlow = 0, xhigh = 1560,
                                   ylow = 0, yhigh = 1,
                                   mods = c("exp", "weibull", "gamma", "lnorm")) 

str(surv_first_flare_treat, max.level = 2) 
names(surv_first_flare_treat$model.objects$models) # List of models fitted 
surv_first_flare_treat$AIC 
surv_first_flare_treat$BIC 

surv_first_flare_control_original <- fit.fun(time = "time", status = "status", data = IPD_first_flare_data_control,
                                               xlow = 0, xhigh = 53,
                                               ylow = 0.75, yhigh = 1,
                                               mods = c("exp", "weibull", "gamma", "lnorm")
) 

surv_first_flare_control <- fit.fun( time = "time", status = "status", data = IPD_first_flare_data_control,
                                       extrapolate = TRUE, times=seq(0, 1560, by = 100),
                                       xlow = 0, xhigh = 1560,
                                       ylow = 0, yhigh = 1,
                                       mods = c("exp", "weibull", "gamma", "lnorm")) 

surv_first_flare_control$AIC 
surv_first_flare_control$BIC 

# Following ERG's suggestions, choose Weibull distribution 
# KM plot
plot(km_treat,
     ylab = "Survival Probability",
     xlab = "Time",
     conf.int = FALSE,
     ylim = c(0, 1),
     xlim = c(0, 480),
     main = "True vs Fitted OS",
     lwd = 2,
     col = "black")
# Add Weibull line
weibull_model_treat <- surv_first_flare_treat$model.objects$models[["Weibull (AFT)"]]
weibull_surv_treat <- summary(weibull_model_treat, type = "survival", t = seq(0, 1560, by = 100), ci = FALSE)
weibull_surv_treat_df <- as.data.frame(weibull_surv_treat[[1]])
lines(weibull_surv_treat_df$time, weibull_surv_treat_df$est, col = 2, lty = 2, lwd = 2)
# Legend
legend("topright", legend = c("Kaplan-Meier", "Weibull Fit"),
       col = c("black", "red"), lty = c(1, 2), lwd = 2)

plot(km_control,
     ylab = "Survival Probability",
     xlab = "Time",
     conf.int = FALSE,
     ylim = c(0, 1),
     xlim = c(0, 480),
     main = "True vs Fitted OS",
     lwd = 2,
     col = "black")
weibull_model_control <- surv_first_flare_control$model.objects$models[["Weibull (AFT)"]]
weibull_surv_control <- summary(weibull_model_control, type = "survival", t = seq(0, 1560, by = 100), ci = FALSE)
weibull_surv_control_df <- as.data.frame(weibull_surv_control[[1]])
lines(weibull_surv_control_df$time, weibull_surv_control_df$est, col = 2, lty = 2, lwd = 2)
legend("topright", legend = c("Kaplan-Meier", "Weibull Fit"),
       col = c("black", "red"), lty = c(1, 2), lwd = 2)

# Calculate survival probabilities.
surv_prob <- function(model, times = NULL, PA = FALSE, rx = 1) {
  if (PA) {
    surv <- model$mat[[rx]][,-1]
  } else {
    surv <- summary(model, t = times, ci = F)[[rx]]$est
  }
  return(surv)
}

surv_probs_first_flare_treat <- surv_prob(weibull_model_treat, times = seq(0, 52*30, by = 52/12))
surv_probs_first_flare_control <- surv_prob(weibull_model_control, times = seq(0, 52*30, by = 52/12))

# Calculate transition probabilities.
trans_prob <- function(surv) {
  # Calculate discrete hazard (i.e., conditional probability of event in each interval)
  d_surv <- surv[-1] / surv[-length(surv)]
  
  # Handle potential division by zero or undefined values
  d_surv[is.na(d_surv)] <- 1  # i.e., assume complete drop (zero survival) → full transition
  
  # Ensure discrete hazards don't exceed 1 (numerical stability)
  d_surv[d_surv > 1] <- 1
  
  # Transition probability is 1 - (S_t+1 / S_t)
  t.p <- 1 - d_surv
  
  # Set negative probabilities (if any) to 0 and issue a message
  if (any(t.p < 0)) {
    message("Negative transition probabilities were set to 0.")
    t.p[t.p < 0] <- 0
  }
  
  return(t.p)
}

trans_probs_first_flare_treat <- trans_prob(surv_probs_first_flare_treat)
trans_probs_first_flare_control <- trans_prob(surv_probs_first_flare_control)

# Display result
df_probs_first_flare_treat <- data.frame(Cycle = 1:length(trans_probs_first_flare_treat),
                                           Trans_Prob = trans_probs_first_flare_treat)

df_probs_first_flare_control <- data.frame(Cycle = 1:length(trans_probs_first_flare_control),
                                             Trans_Prob = trans_probs_first_flare_control)


