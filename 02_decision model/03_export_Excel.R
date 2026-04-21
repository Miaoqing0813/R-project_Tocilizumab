#######################################
# Export results to Excel workbook
#######################################

# Create a new workbook
wb <- createWorkbook()

# Define reusable styles
headerStyle <- createStyle(
  fontSize = 11, fontColour = "#FFFFFF", halign = "center",
  fgFill = "#4F81BD", border = "TopBottomLeftRight",
  borderColour = "#000000", textDecoration = "bold"
)

sectionStyle <- createStyle(
  fontSize = 12, fontColour = "#FFFFFF", halign = "left",
  fgFill = "#366092", border = "TopBottomLeftRight",
  borderColour = "#000000", textDecoration = "bold"
)

titleStyle <- createStyle(
  fontSize = 14, fontColour = "#FFFFFF", halign = "center",
  fgFill = "#1F4E78", border = "TopBottomLeftRight",
  borderColour = "#000000", textDecoration = "bold"
)

costStyle <- createStyle(numFmt = "£#,##0.00")
icerStyle <- createStyle(numFmt = "£#,##0")
qalyStyle <- createStyle(numFmt = "0.0000")
countStyle <- createStyle(numFmt = "0.000")
rateStyle <- createStyle(numFmt = "0.00")
percentStyle <- createStyle(numFmt = "0.00")

# SHEET 1: EXECUTIVE SUMMARY
#######################################

addWorksheet(wb, "Executive Summary")
current_row <- 1

# Title
writeData(wb, "Executive Summary", "PSA RESULTS - EXECUTIVE SUMMARY",
          startRow = current_row, startCol = 1)
addStyle(wb, "Executive Summary", titleStyle, rows = current_row, cols = 1:5, gridExpand = TRUE)
mergeCells(wb, "Executive Summary", rows = current_row, cols = 1:5)
current_row <- current_row + 2

# Key Results
writeData(wb, "Executive Summary", "KEY COST-EFFECTIVENESS RESULTS",
          startRow = current_row, startCol = 1)
addStyle(wb, "Executive Summary", sectionStyle, rows = current_row, cols = 1:5, gridExpand = TRUE)
current_row <- current_row + 1

key_results <- data.frame(
  Outcome = c("Incremental Cost", "Incremental QALYs", "ICER (£/QALY)", 
              "Probability CE at £20,000", "Probability CE at £30,000"),
  Mean = c(
    mean(inc_costs_psa),
    mean(inc_qalys_psa),
    mean(psa_icer[is.finite(psa_icer)]),
    mean((inc_costs_psa / inc_qalys_psa) <= 20000 | (inc_costs_psa <= 0 & inc_qalys_psa >= 0)) * 100,
    mean((inc_costs_psa / inc_qalys_psa) <= 30000 | (inc_costs_psa <= 0 & inc_qalys_psa >= 0)) * 100
  ),
  Lower_95CI = c(
    quantile(inc_costs_psa, 0.025),
    quantile(inc_qalys_psa, 0.025),
    quantile(psa_icer[is.finite(psa_icer)], 0.025),
    NA,
    NA
  ),
  Upper_95CI = c(
    quantile(inc_costs_psa, 0.975),
    quantile(inc_qalys_psa, 0.975),
    quantile(psa_icer[is.finite(psa_icer)], 0.975),
    NA,
    NA
  ),
  stringsAsFactors = FALSE
)

key_results$Mean <- round(key_results$Mean, 2)
key_results$Lower_95CI <- round(key_results$Lower_95CI, 2)
key_results$Upper_95CI <- round(key_results$Upper_95CI, 2)

writeData(wb, "Executive Summary", key_results, startRow = current_row, startCol = 1)
addStyle(wb, "Executive Summary", headerStyle, rows = current_row,
         cols = 1:4, gridExpand = TRUE)

# Format specific rows
addStyle(wb, "Executive Summary", costStyle, rows = current_row + 1, cols = 2:4, gridExpand = TRUE)
addStyle(wb, "Executive Summary", qalyStyle, rows = current_row + 2, cols = 2:4, gridExpand = TRUE)
addStyle(wb, "Executive Summary", icerStyle, rows = current_row + 3, cols = 2:4, gridExpand = TRUE)
addStyle(wb, "Executive Summary", percentStyle, rows = (current_row + 4):(current_row + 5), 
         cols = 2, gridExpand = TRUE)

setColWidths(wb, "Executive Summary", cols = 1:4, widths = "auto")

# SHEET 2: DETAILED PSA RESULTS
#######################################

addWorksheet(wb, "Detailed PSA Results")
current_row <- 1

# Section 1: Summary Statistics
writeData(wb, "Detailed PSA Results", "COST AND QALY SUMMARY STATISTICS",
          startRow = current_row, startCol = 1)
addStyle(wb, "Detailed PSA Results", sectionStyle, rows = current_row, 
         cols = 1:6, gridExpand = TRUE)
current_row <- current_row + 1

summary_stats <- data.frame(
  Outcome = c("Total Cost - TCZ+GC", "Total Cost - GC alone", "Incremental Cost",
              "Total QALYs - TCZ+GC", "Total QALYs - GC alone", "Incremental QALYs",
              "ICER (£/QALY)"),
  Mean = c(
    mean(psa_total_costs[, 1]),
    mean(psa_total_costs[, 2]),
    mean(inc_costs_psa),
    mean(psa_total_qalys[, 1]),
    mean(psa_total_qalys[, 2]),
    mean(inc_qalys_psa),
    mean(psa_icer[is.finite(psa_icer)])
  ),
  Median = c(
    median(psa_total_costs[, 1]),
    median(psa_total_costs[, 2]),
    median(inc_costs_psa),
    median(psa_total_qalys[, 1]),
    median(psa_total_qalys[, 2]),
    median(inc_qalys_psa),
    median(psa_icer[is.finite(psa_icer)])
  ),
  SD = c(
    sd(psa_total_costs[, 1]),
    sd(psa_total_costs[, 2]),
    sd(inc_costs_psa),
    sd(psa_total_qalys[, 1]),
    sd(psa_total_qalys[, 2]),
    sd(inc_qalys_psa),
    sd(psa_icer[is.finite(psa_icer)])
  ),
  Lower_95CI = c(
    quantile(psa_total_costs[, 1], 0.025),
    quantile(psa_total_costs[, 2], 0.025),
    quantile(inc_costs_psa, 0.025),
    quantile(psa_total_qalys[, 1], 0.025),
    quantile(psa_total_qalys[, 2], 0.025),
    quantile(inc_qalys_psa, 0.025),
    quantile(psa_icer[is.finite(psa_icer)], 0.025)
  ),
  Upper_95CI = c(
    quantile(psa_total_costs[, 1], 0.975),
    quantile(psa_total_costs[, 2], 0.975),
    quantile(inc_costs_psa, 0.975),
    quantile(psa_total_qalys[, 1], 0.975),
    quantile(psa_total_qalys[, 2], 0.975),
    quantile(inc_qalys_psa, 0.975),
    quantile(psa_icer[is.finite(psa_icer)], 0.975)
  ),
  stringsAsFactors = FALSE
)

writeData(wb, "Detailed PSA Results", summary_stats, startRow = current_row, startCol = 1)
addStyle(wb, "Detailed PSA Results", headerStyle, rows = current_row,
         cols = 1:6, gridExpand = TRUE)

# Format costs
addStyle(wb, "Detailed PSA Results", costStyle,
         rows = (current_row + 1):(current_row + 3),
         cols = 2:6, gridExpand = TRUE)

# Format QALYs
addStyle(wb, "Detailed PSA Results", qalyStyle,
         rows = (current_row + 4):(current_row + 6),
         cols = 2:6, gridExpand = TRUE)

# Format ICER
addStyle(wb, "Detailed PSA Results", icerStyle,
         rows = current_row + 7,
         cols = 2:6, gridExpand = TRUE)

current_row <- current_row + nrow(summary_stats) + 3

# Section 2: Cost-Effectiveness Probabilities
writeData(wb, "Detailed PSA Results", "COST-EFFECTIVENESS PROBABILITIES",
          startRow = current_row, startCol = 1)
addStyle(wb, "Detailed PSA Results", sectionStyle, rows = current_row, 
         cols = 1:3, gridExpand = TRUE)
current_row <- current_row + 1

wtp_thresholds <- c(0, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 50000)
prob_ce_table <- data.frame(
  WTP_Threshold = paste0("£", format(wtp_thresholds, big.mark = ",")),
  Probability = numeric(length(wtp_thresholds)),
  Percentage = numeric(length(wtp_thresholds)),
  stringsAsFactors = FALSE
)

for (i in 1:length(wtp_thresholds)) {
  wtp <- wtp_thresholds[i]
  if (wtp == 0) {
    prob <- mean(inc_costs_psa <= 0 & inc_qalys_psa >= 0)
  } else {
    prob <- mean((inc_costs_psa / inc_qalys_psa) <= wtp | (inc_costs_psa <= 0 & inc_qalys_psa >= 0))
  }
  prob_ce_table$Probability[i] <- round(prob, 4)
  prob_ce_table$Percentage[i] <- round(prob * 100, 2)
}

writeData(wb, "Detailed PSA Results", prob_ce_table, startRow = current_row, startCol = 1)
addStyle(wb, "Detailed PSA Results", headerStyle, rows = current_row,
         cols = 1:3, gridExpand = TRUE)

addStyle(wb, "Detailed PSA Results", qalyStyle,
         rows = (current_row + 1):(current_row + nrow(prob_ce_table)),
         cols = 2, gridExpand = TRUE)

addStyle(wb, "Detailed PSA Results", percentStyle,
         rows = (current_row + 1):(current_row + nrow(prob_ce_table)),
         cols = 3, gridExpand = TRUE)

setColWidths(wb, "Detailed PSA Results", cols = 1:6, widths = "auto")

# SHEET 3: ADVERSE EVENTS - INCIDENTS
#######################################

addWorksheet(wb, "AE Incidents")
current_row <- 1

n_events <- dim(psa_incident_events)[2]
event_names <- dimnames(psa_incident_events)$event

writeData(wb, "AE Incidents", "TOTAL INCIDENT EVENTS",
          startRow = current_row, startCol = 1)
addStyle(wb, "AE Incidents", sectionStyle, rows = current_row, 
         cols = 1:8, gridExpand = TRUE)
current_row <- current_row + 1

ae_incidents <- data.frame(
  Event = event_names,
  TCZ_GC_Mean = numeric(n_events),
  TCZ_GC_LCI = numeric(n_events),
  TCZ_GC_UCI = numeric(n_events),
  GC_Alone_Mean = numeric(n_events),
  GC_Alone_LCI = numeric(n_events),
  GC_Alone_UCI = numeric(n_events),
  Difference = numeric(n_events),
  stringsAsFactors = FALSE
)

for (i in 1:n_events) {
  ae_incidents$TCZ_GC_Mean[i] <- mean(psa_incident_events[, i, 1])
  ae_incidents$TCZ_GC_LCI[i] <- quantile(psa_incident_events[, i, 1], 0.025)
  ae_incidents$TCZ_GC_UCI[i] <- quantile(psa_incident_events[, i, 1], 0.975)
  ae_incidents$GC_Alone_Mean[i] <- mean(psa_incident_events[, i, 2])
  ae_incidents$GC_Alone_LCI[i] <- quantile(psa_incident_events[, i, 2], 0.025)
  ae_incidents$GC_Alone_UCI[i] <- quantile(psa_incident_events[, i, 2], 0.975)
  ae_incidents$Difference[i] <- ae_incidents$TCZ_GC_Mean[i] - ae_incidents$GC_Alone_Mean[i]
}

ae_incidents[, 2:8] <- round(ae_incidents[, 2:8], 3)

writeData(wb, "AE Incidents", ae_incidents, startRow = current_row, startCol = 1)
addStyle(wb, "AE Incidents", headerStyle, rows = current_row,
         cols = 1:8, gridExpand = TRUE)

addStyle(wb, "AE Incidents", countStyle,
         rows = (current_row + 1):(current_row + n_events),
         cols = 2:8, gridExpand = TRUE)

setColWidths(wb, "AE Incidents", cols = 1:8, widths = "auto")

# SHEET 4: ADVERSE EVENTS - RATES PER 100 PY
#######################################

addWorksheet(wb, "AE Rates per 100 PY")
current_row <- 1

writeData(wb, "AE Rates per 100 PY", "EVENT RATES PER 100 PATIENT-YEARS",
          startRow = current_row, startCol = 1)
addStyle(wb, "AE Rates per 100 PY", sectionStyle, rows = current_row, 
         cols = 1:8, gridExpand = TRUE)
current_row <- current_row + 1

ae_rates <- data.frame(
  Event = event_names,
  TCZ_GC_Mean = numeric(n_events),
  TCZ_GC_LCI = numeric(n_events),
  TCZ_GC_UCI = numeric(n_events),
  GC_Alone_Mean = numeric(n_events),
  GC_Alone_LCI = numeric(n_events),
  GC_Alone_UCI = numeric(n_events),
  Difference = numeric(n_events),
  stringsAsFactors = FALSE
)

for (i in 1:n_events) {
  ae_rates$TCZ_GC_Mean[i] <- mean(psa_event_rates_per_100py[, i, 1])
  ae_rates$TCZ_GC_LCI[i] <- quantile(psa_event_rates_per_100py[, i, 1], 0.025)
  ae_rates$TCZ_GC_UCI[i] <- quantile(psa_event_rates_per_100py[, i, 1], 0.975)
  ae_rates$GC_Alone_Mean[i] <- mean(psa_event_rates_per_100py[, i, 2])
  ae_rates$GC_Alone_LCI[i] <- quantile(psa_event_rates_per_100py[, i, 2], 0.025)
  ae_rates$GC_Alone_UCI[i] <- quantile(psa_event_rates_per_100py[, i, 2], 0.975)
  ae_rates$Difference[i] <- ae_rates$TCZ_GC_Mean[i] - ae_rates$GC_Alone_Mean[i]
}

ae_rates[, 2:8] <- round(ae_rates[, 2:8], 2)

writeData(wb, "AE Rates per 100 PY", ae_rates, startRow = current_row, startCol = 1)
addStyle(wb, "AE Rates per 100 PY", headerStyle, rows = current_row,
         cols = 1:8, gridExpand = TRUE)

addStyle(wb, "AE Rates per 100 PY", rateStyle,
         rows = (current_row + 1):(current_row + n_events),
         cols = 2:8, gridExpand = TRUE)

# SHEET 5: PER-EVENT COSTS
#######################################

addWorksheet(wb, "Per-Event Costs")
current_row <- 1

writeData(wb, "Per-Event Costs", "LIFETIME PER-EVENT COSTS",
          startRow = current_row, startCol = 1)
addStyle(wb, "Per-Event Costs", sectionStyle, rows = current_row, 
         cols = 1:8, gridExpand = TRUE)
current_row <- current_row + 1

per_event_costs_summary <- data.frame(
  Event = event_names,
  TCZ_GC_Mean = numeric(n_events),
  TCZ_GC_LCI = numeric(n_events),
  TCZ_GC_UCI = numeric(n_events),
  GC_Alone_Mean = numeric(n_events),
  GC_Alone_LCI = numeric(n_events),
  GC_Alone_UCI = numeric(n_events),
  Difference = numeric(n_events),
  stringsAsFactors = FALSE
)

for (i in 1:n_events) {
  per_event_costs_summary$TCZ_GC_Mean[i] <- mean(psa_per_event_costs[, i, 1])
  per_event_costs_summary$TCZ_GC_LCI[i] <- quantile(psa_per_event_costs[, i, 1], 0.025)
  per_event_costs_summary$TCZ_GC_UCI[i] <- quantile(psa_per_event_costs[, i, 1], 0.975)
  per_event_costs_summary$GC_Alone_Mean[i] <- mean(psa_per_event_costs[, i, 2])
  per_event_costs_summary$GC_Alone_LCI[i] <- quantile(psa_per_event_costs[, i, 2], 0.025)
  per_event_costs_summary$GC_Alone_UCI[i] <- quantile(psa_per_event_costs[, i, 2], 0.975)
  per_event_costs_summary$Difference[i] <- per_event_costs_summary$TCZ_GC_Mean[i] - 
    per_event_costs_summary$GC_Alone_Mean[i]
}

per_event_costs_summary[, 2:8] <- round(per_event_costs_summary[, 2:8], 2)

writeData(wb, "Per-Event Costs", per_event_costs_summary, startRow = current_row, startCol = 1)
addStyle(wb, "Per-Event Costs", headerStyle, rows = current_row,
         cols = 1:8, gridExpand = TRUE)

addStyle(wb, "Per-Event Costs", costStyle,
         rows = (current_row + 1):(current_row + n_events),
         cols = 2:8, gridExpand = TRUE)

setColWidths(wb, "Per-Event Costs", cols = 1:8, widths = "auto")

# SHEET 6: PER-EVENT QALY LOSSES
#######################################

addWorksheet(wb, "Per-Event QALY Losses")
current_row <- 1

writeData(wb, "Per-Event QALY Losses", "LIFETIME PER-EVENT QALY LOSSES",
          startRow = current_row, startCol = 1)
addStyle(wb, "Per-Event QALY Losses", sectionStyle, rows = current_row, 
         cols = 1:8, gridExpand = TRUE)
current_row <- current_row + 1

per_event_qalys_summary <- data.frame(
  Event = event_names,
  TCZ_GC_Mean = numeric(n_events),
  TCZ_GC_LCI = numeric(n_events),
  TCZ_GC_UCI = numeric(n_events),
  GC_Alone_Mean = numeric(n_events),
  GC_Alone_LCI = numeric(n_events),
  GC_Alone_UCI = numeric(n_events),
  Difference = numeric(n_events),
  stringsAsFactors = FALSE
)

for (i in 1:n_events) {
  per_event_qalys_summary$TCZ_GC_Mean[i] <- mean(psa_per_event_qalys[, i, 1])
  per_event_qalys_summary$TCZ_GC_LCI[i] <- quantile(psa_per_event_qalys[, i, 1], 0.025)
  per_event_qalys_summary$TCZ_GC_UCI[i] <- quantile(psa_per_event_qalys[, i, 1], 0.975)
  per_event_qalys_summary$GC_Alone_Mean[i] <- mean(psa_per_event_qalys[, i, 2])
  per_event_qalys_summary$GC_Alone_LCI[i] <- quantile(psa_per_event_qalys[, i, 2], 0.025)
  per_event_qalys_summary$GC_Alone_UCI[i] <- quantile(psa_per_event_qalys[, i, 2], 0.975)
  per_event_qalys_summary$Difference[i] <- per_event_qalys_summary$TCZ_GC_Mean[i] - 
    per_event_qalys_summary$GC_Alone_Mean[i]
}

per_event_qalys_summary[, 2:8] <- round(per_event_qalys_summary[, 2:8], 4)

writeData(wb, "Per-Event QALY Losses", per_event_qalys_summary, startRow = current_row, startCol = 1)
addStyle(wb, "Per-Event QALY Losses", headerStyle, rows = current_row,
         cols = 1:8, gridExpand = TRUE)

addStyle(wb, "Per-Event QALY Losses", qalyStyle,
         rows = (current_row + 1):(current_row + n_events),
         cols = 2:8, gridExpand = TRUE)

setColWidths(wb, "Per-Event QALY Losses", cols = 1:8, widths = "auto")

# SHEET 7: RAW SIMULATION DATA
#######################################

addWorksheet(wb, "Raw Simulation Data")

psa_raw_data <- data.frame(
  Simulation = 1:n_sim,
  Cost_TCZ_GC = round(psa_total_costs[, 1], 2),
  Cost_GC_Alone = round(psa_total_costs[, 2], 2),
  Inc_Cost = round(inc_costs_psa, 2),
  QALY_TCZ_GC = round(psa_total_qalys[, 1], 4),
  QALY_GC_Alone = round(psa_total_qalys[, 2], 4),
  Inc_QALY = round(inc_qalys_psa, 4),
  ICER = round(psa_icer, 2)
)

writeData(wb, "Raw Simulation Data", psa_raw_data, startRow = 1, startCol = 1)
addStyle(wb, "Raw Simulation Data", headerStyle, rows = 1, cols = 1:10, gridExpand = TRUE)

# Format columns
addStyle(wb, "Raw Simulation Data", costStyle, 
         rows = 2:(n_sim + 1), cols = c(2, 3, 4), gridExpand = TRUE)
addStyle(wb, "Raw Simulation Data", qalyStyle, 
         rows = 2:(n_sim + 1), cols = c(5, 6, 7), gridExpand = TRUE)
addStyle(wb, "Raw Simulation Data", icerStyle, 
         rows = 2:(n_sim + 1), cols = 8, gridExpand = TRUE)
addStyle(wb, "Raw Simulation Data", rateStyle, 
         rows = 2:(n_sim + 1), cols = c(9, 10), gridExpand = TRUE)

setColWidths(wb, "Raw Simulation Data", cols = 1:10, widths = "auto")

# SHEET 8: CEAC DATA
#######################################

addWorksheet(wb, "CEAC Data")

wtp_range <- seq(0, 50000, by = 1000)
prob_ce_ceac <- numeric(length(wtp_range))

for (i in 1:length(wtp_range)) {
  wtp <- wtp_range[i]
  if (wtp == 0) {
    prob_ce_ceac[i] <- mean(inc_costs_psa <= 0 & inc_qalys_psa >= 0)
  } else {
    prob_ce_ceac[i] <- mean((inc_costs_psa / inc_qalys_psa) <= wtp | (inc_costs_psa <= 0 & inc_qalys_psa >= 0))
  }
}

ceac_data <- data.frame(
  WTP_Threshold = wtp_range,
  Probability = round(prob_ce_ceac, 4),
  Percentage = round(prob_ce_ceac * 100, 2)
)

writeData(wb, "CEAC Data", ceac_data, startRow = 1, startCol = 1)
addStyle(wb, "CEAC Data", headerStyle, rows = 1, cols = 1:3, gridExpand = TRUE)

addStyle(wb, "CEAC Data", icerStyle, rows = 2:(nrow(ceac_data) + 1), 
         cols = 1, gridExpand = TRUE)
addStyle(wb, "CEAC Data", qalyStyle, rows = 2:(nrow(ceac_data) + 1), 
         cols = 2, gridExpand = TRUE)
addStyle(wb, "CEAC Data", percentStyle, rows = 2:(nrow(ceac_data) + 1), 
         cols = 3, gridExpand = TRUE)

setColWidths(wb, "CEAC Data", cols = 1:3, widths = "auto")

# SHEET 9: CE PLANE DATA
#######################################

addWorksheet(wb, "CE Plane Data")

ce_plane_data <- data.frame(
  Simulation = 1:n_sim,
  Inc_QALYs = round(inc_qalys_psa, 4),
  Inc_Costs = round(inc_costs_psa, 2),
  ICER = round(psa_icer, 2),
  Quadrant = character(n_sim),
  stringsAsFactors = FALSE
)

ce_plane_data$Quadrant <- ifelse(
  inc_qalys_psa >= 0 & inc_costs_psa >= 0, "NE - More effective, more costly",
  ifelse(inc_qalys_psa >= 0 & inc_costs_psa < 0, "SE - Dominant",
         ifelse(inc_qalys_psa < 0 & inc_costs_psa >= 0, "NW - Dominated",
                "SW - Less effective, less costly")))

writeData(wb, "CE Plane Data", ce_plane_data, startRow = 1, startCol = 1)
addStyle(wb, "CE Plane Data", headerStyle, rows = 1, cols = 1:5, gridExpand = TRUE)

addStyle(wb, "CE Plane Data", qalyStyle, rows = 2:(n_sim + 1), cols = 2, gridExpand = TRUE)
addStyle(wb, "CE Plane Data", costStyle, rows = 2:(n_sim + 1), cols = 3, gridExpand = TRUE)
addStyle(wb, "CE Plane Data", icerStyle, rows = 2:(n_sim + 1), cols = 4, gridExpand = TRUE)

setColWidths(wb, "CE Plane Data", cols = 1:5, widths = "auto")

# SAVE WORKBOOK
#######################################

output_file <- "03_results/PSA_results.xlsx"
saveWorkbook(wb, output_file, overwrite = TRUE)

cat("\n==================================================\n")
cat(paste("PSA results successfully saved to:", output_file, "\n"))
cat("==================================================\n\n")
cat("Workbook contains the following sheets:\n")
cat("  1. Executive Summary - Key results at a glance\n")
cat("  2. Detailed PSA Results - Full statistical summary\n")
cat("  3. AE Incidents - Total incident events with CIs\n")
cat("  4. AE Rates per 100 PY - Event rates per 100 patient-years\n")
cat("  5. Per-Event Costs - Lifetime costs by event type\n")
cat("  6. Per-Event QALY Losses - QALY impacts by event type\n")
cat("  7. Raw Simulation Data - Individual simulation results\n")
cat("  8. CEAC Data - Cost-effectiveness acceptability curve\n")
cat("  9. CE Plane Data - Cost-effectiveness plane with quadrants\n\n")

