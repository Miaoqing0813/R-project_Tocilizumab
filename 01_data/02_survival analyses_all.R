##############################################
# Transition probabilities to first flare #
##############################################

#' Generate transition probabilities for first flare
#'
#' @param subgroup Character: "itt", "newlydiag", or "relapsed"
#' @param treatment_duration_weeks Numeric: Treatment duration in weeks (e.g., 52 or 104)
#' @param time_horizon_weeks Numeric: Total time horizon for extrapolation in weeks (default: 5*52)
#' @param time_horizon_months Numeric: Number of monthly cycles (default: 5*12)
#' @param data_path Character: Path to the directory containing KM curve images
#' @param generate_plots Logical: Whether to generate diagnostic plots (default: FALSE)
#' @param verbose Logical: Whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#'   - trans_probs_treat: Data frame with treatment arm transition probabilities
#'   - trans_probs_control: Data frame with control arm transition probabilities
#'   - models: List of fitted survival models (within and beyond treatment duration)
#'   - ipd_data: List of IPD datasets used for fitting

generate_first_flare_transition_probs <- function(
    subgroup = c("itt", "newlydiag", "relapsed"),
    treatment_duration_weeks = 52,
    time_horizon_weeks = 5 * 52,
    time_horizon_months = 5 * 12,
    dist_beyond = "exp", 
    data_path = here::here(getwd(), "01_data/Raw data"),
    generate_plots = FALSE,
    verbose = TRUE
) {
  
  # Validate inputs
  subgroup <- match.arg(subgroup)
  
  if (verbose) {
    cat("\n=== Starting Survival Analysis ===\n")
    cat("Subgroup:", subgroup, "\n")
    cat("Treatment duration:", treatment_duration_weeks, "weeks\n")
    cat("Time horizon:", time_horizon_weeks/52, "years\n")
    cat("Number of cycles:", time_horizon_months, "\n\n")
  }
  
  # Define file paths based on subgroup
  km_file_map <- list(
    itt = "Time to first flare_ITT.png",
    newlydiag = "Time to first flare_newly diagnosed.png",
    relapsed = "Time to first flare_relapsed.png"
  )
  
  km_file <- file.path(data_path, km_file_map[[subgroup]])
  
  if (!file.exists(km_file)) {
    stop("KM curve file not found: ", km_file)
  }
  
  # Define curve indices and at-risk numbers based on subgroup
  curve_config <- list(
    itt = list(
      treat_curve = 3,
      control_curve = 1,
      n_risk_treat = c(100, 93, 88, 85, 85, 81, 77, 74, 71, 69, 67, 64, 63, 5),
      n_risk_control = c(51, 48, 44, 41, 38, 35, 32, 30, 28, 25, 22, 17, 15, 0)
    ),
    newlydiag = list(
      treat_curve = 1,
      control_curve = 4,
      n_risk_treat = c(47, 43, 42, 40, 40, 40, 38, 36, 35, 34, 33, 33, 32, 4),
      n_risk_control = c(23, 21, 21, 19, 18, 17, 16, 14, 14, 13, 13, 11, 10, 0)
    ),
    relapsed = list(
      treat_curve = 3,
      control_curve = 4,
      n_risk_treat = c(53, 50, 46, 45, 45, 41, 39, 38, 36, 35, 34, 31, 31, 1),
      n_risk_control = c(28, 27, 23, 22, 20, 18, 16, 16, 14, 12, 9, 6, 5, 0)
    )
  )
  
  config <- curve_config[[subgroup]]
  
  # Step 1: Digitize KM curves
  if (verbose) cat("Step 1: Digitizing KM curves...\n")
  
  first_flare_data <- SurvdigitizeR::survival_digitize(
    img_path = km_file,
    num_curves = 4,
    x_start = 0,
    x_end = 52,
    x_increment = 4,
    y_start = 0,
    y_end = 100,
    y_increment = 20,
    y_text_vertical = TRUE
  )
  
  # Extract treatment and control curves
  first_flare_data_treat <- first_flare_data[which(first_flare_data$curve == config$treat_curve), c(2, 3)]
  first_flare_data_control <- first_flare_data[which(first_flare_data$curve == config$control_curve), c(2, 3)]
  
  # Step 2: Preprocess data
  if (verbose) cat("Step 2: Preprocessing data...\n")
  
  weeks <- c(0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52)
  
  pre_first_flare_data_treat <- IPDfromKM::preprocess(
    first_flare_data_treat,
    trisk = weeks,
    nrisk = config$n_risk_treat,
    totalpts = NULL,
    maxy = 100
  )
  
  pre_first_flare_data_control <- IPDfromKM::preprocess(
    first_flare_data_control,
    trisk = weeks,
    nrisk = config$n_risk_control,
    totalpts = NULL,
    maxy = 100
  )
  
  # Step 3: Generate IPD
  if (verbose) cat("Step 3: Generating individual patient data...\n")
  
  IPD_first_flare_data_treat <- IPDfromKM::getIPD(prep = pre_first_flare_data_treat, armID = 1)$IPD
  IPD_first_flare_data_control <- IPDfromKM::getIPD(prep = pre_first_flare_data_control, armID = 0)$IPD
  
  # Step 4: Fit survival models - WITHIN treatment period
  if (verbose) cat("Step 4: Fitting survival models (within treatment period)...\n")
  
  # Fit multiple survival models with different distributions
  survival_models_treat <- list(
    exp = flexsurv::flexsurvreg(survival::Surv(time, status) ~ 1, 
                                data = IPD_first_flare_data_treat, 
                                dist = "exp"),
    weibull = flexsurv::flexsurvreg(survival::Surv(time, status) ~ 1, 
                                    data = IPD_first_flare_data_treat, 
                                    dist = "weibull"), 
    gamma = flexsurv::flexsurvreg(survival::Surv(time, status) ~ 1, 
                                  data = IPD_first_flare_data_treat, 
                                  dist = "gamma"),
    lnorm = flexsurv::flexsurvreg(survival::Surv(time, status) ~ 1, 
                                  data = IPD_first_flare_data_treat, 
                                  dist = "lnorm")
  )
  
  survival_models_control <- list(
    exp = flexsurv::flexsurvreg(survival::Surv(time, status) ~ 1, 
                                data = IPD_first_flare_data_control, 
                                dist = "exp"),
    weibull = flexsurv::flexsurvreg(survival::Surv(time, status) ~ 1, 
                                    data = IPD_first_flare_data_control, 
                                    dist = "weibull"),
    gamma = flexsurv::flexsurvreg(survival::Surv(time, status) ~ 1, 
                                  data = IPD_first_flare_data_control, 
                                  dist = "gamma"),
    lnorm = flexsurv::flexsurvreg(survival::Surv(time, status) ~ 1, 
                                  data = IPD_first_flare_data_control, 
                                  dist = "lnorm")
  )
  
  # Select best-fitting models based on subgroup during the treatment period (based on original HTA choices)
  model_selection <- list(
    itt = list(treat_within = "weibull", control_within = "exp"),
    newlydiag = list(treat_within = "exp", control_within = "exp"),
    relapsed = list(treat_within = "exp", control_within = "lnorm")
  )
  
  selected_models <- model_selection[[subgroup]]
  
  model_treat_within <- survival_models_treat[[selected_models$treat_within]]
  model_control_within <- survival_models_control[[selected_models$control_within]]
  
  if (verbose) {
    cat("  Treatment arm (within): Using", selected_models$treat_within, "\n")
    cat("  Control arm (within): Using", selected_models$control_within, "\n")
  }
  
  # Step 5: Fit survival models - BEYOND treatment period
  if (verbose) cat("Step 5: Fitting survival models (beyond treatment period)...\n")
  
  # Choose distributions for extrapolation 
  model_treat_beyond <- survival_models_treat[[dist_beyond]]
  model_control_beyond <- survival_models_control[[dist_beyond]]
  
  if (verbose) {
    cat("  Both arms (beyond): Using", dist_beyond, "\n")
  }
  
  # Step 6: Calculate transition probabilities
  if (verbose) cat("Step 6: Calculating transition probabilities...\n")
  
  trans_probs_treat <- calculate_transition_probs(
    short_term_model = model_treat_within,
    long_term_model = model_treat_beyond,
    treatment_duration_weeks = treatment_duration_weeks,
    time_horizon_weeks = time_horizon_weeks,
    time_horizon_months = time_horizon_months
  )
  
  trans_probs_control <- calculate_transition_probs(
    short_term_model = model_control_within,
    long_term_model = model_control_beyond,
    treatment_duration_weeks = treatment_duration_weeks,
    time_horizon_weeks = time_horizon_weeks,
    time_horizon_months = time_horizon_months
  )
  
  # Step 7: Generate plots if requested
  if (generate_plots) {
    if (verbose) cat("Step 7: Generating diagnostic plots...\n")
    
    # KM plots
    km_treat_plot <- survfit2(survival::Surv(time = time, event = status) ~ 1, 
                   data = IPD_first_flare_data_treat) %>%
      ggsurvfit::ggsurvfit() +
      labs(
        x = "Weeks",
        y = paste0("Patients without flare -- Treatment (", subgroup, ")"),
        title = paste0("Kaplan-Meier: Treatment Arm (", subgroup, ")")
      ) +
      ggsurvfit::add_confidence_interval() +
      ggsurvfit::add_risktable()
    
    km_control_plot <- survfit2(survival::Surv(time = time, event = status) ~ 1, 
                   data = IPD_first_flare_data_control) %>%
      ggsurvfit::ggsurvfit() +
      labs(
        x = "Weeks",
        y = paste0("Patients without flare -- Control (", subgroup, ")"),
        title = paste0("Kaplan-Meier: Control Arm (", subgroup, ")")
      ) +
      ggsurvfit::add_confidence_interval() +
      ggsurvfit::add_risktable()
    
    # Fitted survival curves plots
    surv_first_flare_treat_beyond_plot <- fit.fun(
      time = "time",
      status = "status",
      data = IPD_first_flare_data_treat,
      extrapolate = TRUE,
      times = seq(0, time_horizon_weeks, by = 52), 
      xlow = 0, xhigh = time_horizon_weeks,
      ylow = 0, yhigh = 1,
      mods = c("exp", "weibull")
    )
    
    surv_first_flare_control_beyond_plot <- fit.fun(
      time = "time",
      status = "status",
      data = IPD_first_flare_data_control,
      extrapolate = TRUE,
      times = seq(0, time_horizon_weeks, by = 52), 
      xlow = 0, xhigh = time_horizon_weeks,
      ylow = 0, yhigh = 1,
      mods = c("exp", "weibull")
    )
    
    # Transition probability plots
    trans_prob_treat_plot <- plot_transition_probs(
      trans_probs_treat,
      title = paste0("Monthly Transition Probabilities: Treatment (", subgroup, ")")
    )
    
    trans_prob_control_plot <- plot_transition_probs(
      trans_probs_control,
      title = paste0("Monthly Transition Probabilities: Control (", subgroup, ")")
    )
    
  }
  
  if (verbose) cat("\n=== Analysis Complete ===\n\n")
  
  # Return results
  return(list(
    trans_probs_treat = trans_probs_treat,
    trans_probs_control = trans_probs_control,
    models = list(
      treat_within = model_treat_within,
      treat_beyond = model_treat_beyond,
      control_within = model_control_within,
      control_beyond = model_control_beyond
    ),
    ipd_data = list(
      treat = IPD_first_flare_data_treat,
      control = IPD_first_flare_data_control
    ),
    config = list(
      subgroup = subgroup,
      treatment_duration_weeks = treatment_duration_weeks,
      time_horizon_weeks = time_horizon_weeks,
      time_horizon_months = time_horizon_months
    ),
    km_treat_plot = km_treat_plot, 
    km_control_plot = km_control_plot, 
    surv_first_flare_treat_beyond_plot = surv_first_flare_treat_beyond_plot,
    surv_first_flare_control_beyond_plot = surv_first_flare_control_beyond_plot,
    trans_prob_treat_plot = trans_prob_treat_plot, 
    trans_prob_control_plot = trans_prob_control_plot
  ))
}


#' Calculate monthly transition probabilities from survival models
#'
#' @param short_term_model Fitted flexsurvreg model for within-treatment period
#' @param long_term_model Fitted flexsurvreg model for beyond-treatment period
#' @param treatment_duration_weeks Treatment duration in weeks
#' @param time_horizon_weeks Total time horizon in weeks
#' @param time_horizon_months Number of monthly cycles
#'
#' @return Data frame with monthly transition probabilities and confidence intervals

calculate_transition_probs <- function(short_term_model, long_term_model, 
                                       treatment_duration_weeks,
                                       time_horizon_weeks, time_horizon_months) {
  
  # Short-term monthly probabilities (0 to treatment duration)
  monthly_timepoints <- seq(0, treatment_duration_weeks, by = 52/12)
  
  # Get short-term survival probabilities
  short_surv_summary <- summary(short_term_model, type = "survival", t = monthly_timepoints, ci = TRUE)[[1]]
  short_surv_est <- short_surv_summary$est
  short_surv_lower <- short_surv_summary$lcl
  short_surv_upper <- short_surv_summary$ucl
  
  # Calculate monthly transition probabilities for first year
  n_short_months <- length(short_surv_est) - 1
  
  short_trans_prob_est <- 1 - (short_surv_est[-1] / short_surv_est[-length(short_surv_est)])
  short_trans_prob_lower <- 1 - (short_surv_upper[-1] / short_surv_upper[-length(short_surv_upper)])
  short_trans_prob_upper <- 1 - (short_surv_lower[-1] / short_surv_lower[-length(short_surv_lower)])
  
  # Long-term yearly probabilities converted to monthly (beyond treatment duration)
  yearly_timepoints <- seq(treatment_duration_weeks + 52, time_horizon_weeks, by = 52)
  
  # Get long-term survival probabilities
  long_surv_summary <- summary(long_term_model, type = "survival", t = yearly_timepoints, ci = TRUE)[[1]]
  long_surv_est <- pmax(long_surv_summary$est, 0.001)
  long_surv_lower <- pmax(long_surv_summary$lcl, 0.001)
  long_surv_upper <- pmax(long_surv_summary$ucl, 0.001)
  
  # Add baseline from end of short-term period
  endshort_long_surv_est <- c(short_surv_est[length(short_surv_est)], long_surv_est)
  endshort_long_surv_lower <- c(short_surv_lower[length(short_surv_lower)], long_surv_lower)
  endshort_long_surv_upper <- c(short_surv_upper[length(short_surv_upper)], long_surv_upper)
  
  # Calculate yearly transition probabilities for long-term
  n_long_years <- length(endshort_long_surv_est) - 1
  
  long_yearly_trans_est <- 1 - (endshort_long_surv_est[-1] / endshort_long_surv_est[-length(endshort_long_surv_est)])
  long_yearly_trans_lower <- 1 - (endshort_long_surv_lower[-1] / endshort_long_surv_lower[-length(endshort_long_surv_lower)])
  long_yearly_trans_upper <- 1 - (endshort_long_surv_upper[-1] / endshort_long_surv_upper[-length(endshort_long_surv_upper)])
  
  # Ensure transition probabilities are valid (0 to 1)
  long_yearly_trans_est <- pmax(0, pmin(1, long_yearly_trans_est))
  long_yearly_trans_lower <- pmax(0, pmin(1, long_yearly_trans_lower))
  long_yearly_trans_upper <- pmax(0, pmin(1, long_yearly_trans_upper))
  
  # Convert yearly to monthly equivalents
  long_monthly_trans_est <- 1 - (1 - long_yearly_trans_est)^(1/12)
  long_monthly_trans_lower <- 1 - (1 - long_yearly_trans_lower)^(1/12)
  long_monthly_trans_upper <- 1 - (1 - long_yearly_trans_upper)^(1/12)
  
  # Combine the short-term and long-term monthly sequence
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
    start_month <- n_short_months + (i-1)*12 + 1
    end_month <- min(n_short_months + i*12, time_horizon_months)
    
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
    period = c(rep("Short-term", n_short_months), 
               rep("Long-term", time_horizon_months - n_short_months)),
    stringsAsFactors = FALSE
  )
  
  return(trans_prob)
}


#' Plot transition probabilities over time
#'
#' @param trans_prob_df Data frame from calculate_transition_probs
#' @param title Plot title
#' @param max_months Maximum months to display (default: 60)
#'
#' @return ggplot object

plot_transition_probs <- function(trans_prob_df, title="Transition probabilities", max_months = time_horizon_months) {
  df_subset <- trans_prob_df[1:min(max_months, nrow(trans_prob_df)), ]
  
  ggplot2::ggplot(df_subset, ggplot2::aes(x = month, y = trans_prob)) +
    ggplot2::geom_line(color = "blue", size = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = trans_prob_lower, ymax = trans_prob_upper),
                         alpha = 0.3, fill = "blue") +
    ggplot2::labs(
      title = title,
      x = "Month",
      y = "Monthly Transition Probability"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_continuous(breaks = seq(0, max_months, 12)) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12))
}


#' @examples
#' # For newly diagnosed patients with 52-week treatment
#' results_newlydiag_52 <- generate_first_flare_transition_probs(
#'   subgroup = "newlydiag",
#'   treatment_duration_weeks = 52
#' )
#'
#' # For relapsed patients with 104-week treatment
#' results_relapsed_104 <- generate_first_flare_transition_probs(
#'   subgroup = "relapsed",
#'   treatment_duration_weeks = 104
#' )

