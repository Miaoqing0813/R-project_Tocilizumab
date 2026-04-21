####################################################
# EXPORT ALL RESULTS TO EXCEL
####################################################

# Create a new workbook
wb <- createWorkbook()

# Add a worksheet for the cost-effectiveness summary
######################################################
addWorksheet(wb, "CE Summary")
writeData(wb, "CE Summary", summary_table, startRow = 1, startCol = 1)

# Format the CE Summary sheet
headerStyle <- createStyle(fontSize = 11, fontColour = "#FFFFFF", 
                           halign = "center", fgFill = "#4F81BD",
                           border = "TopBottomLeftRight", borderColour = "#000000",
                           textDecoration = "bold")
setColWidths(wb, "CE Summary", cols = 1:ncol(summary_table), widths = "auto")
addStyle(wb, "CE Summary", headerStyle, rows = 1, cols = 1:ncol(summary_table), gridExpand = TRUE)

# Add number formatting for costs
costStyle <- createStyle(numFmt = "£#,##0.00")
addStyle(wb, "CE Summary", costStyle, rows = 2:(nrow(summary_table) + 1), 
         cols = c(4, 5, 6), gridExpand = TRUE)

# Add number formatting for QALYs
qalyStyle <- createStyle(numFmt = "0.0000")
addStyle(wb, "CE Summary", qalyStyle, rows = 2:(nrow(summary_table) + 1), 
         cols = c(7, 8, 9), gridExpand = TRUE)

# Add number formatting for ICER
icerStyle <- createStyle(numFmt = "£#,##0")
addStyle(wb, "CE Summary", icerStyle, rows = 2:(nrow(summary_table) + 1), 
         cols = 10, gridExpand = TRUE)

# Add adverse event tables for each scenarios
###############################################

for (scenario_name in names(scenario_results)) {
  res <- scenario_results[[scenario_name]]
  ae_summary <- res$ae_summary
  
  # Create a safe sheet name (Excel has 31 character limit)
  sheet_name <- paste0("AE_", res$scenario_info$subgroup, "_", 
                       res$scenario_info$treatment_duration_months, "m")
  
  # Truncate if needed
  if (nchar(sheet_name) > 31) {
    sheet_name <- substr(sheet_name, 1, 31)
  }
  
  cat(paste("Adding sheet:", sheet_name, "\n"))
  
  # Add worksheet
  addWorksheet(wb, sheet_name)
  
  current_row <- 1
  
  # Add scenario information header
  writeData(wb, sheet_name, scenario_name, startRow = current_row, startCol = 1)
  titleStyle <- createStyle(fontSize = 14, textDecoration = "bold")
  addStyle(wb, sheet_name, titleStyle, rows = current_row, cols = 1)
  current_row <- current_row + 2
  
  # Add scenario details
  scenario_details <- data.frame(
    Parameter = c("Patient Subgroup", "Treatment Duration"),
    Value = c(res$scenario_info$subgroup, 
              paste(res$scenario_info$treatment_duration_months, "months"))
  )
  writeData(wb, sheet_name, scenario_details, startRow = current_row, startCol = 1)
  addStyle(wb, sheet_name, headerStyle, rows = current_row, cols = 1:2, gridExpand = TRUE)
  current_row <- current_row + nrow(scenario_details) + 2
  
  # Section 1: Total Incident Events
  writeData(wb, sheet_name, "TOTAL INCIDENT EVENTS", startRow = current_row, startCol = 1)
  sectionStyle <- createStyle(fontSize = 12, textDecoration = "bold", fgFill = "#D9E1F2")
  addStyle(wb, sheet_name, sectionStyle, rows = current_row, cols = 1:4, gridExpand = TRUE)
  current_row <- current_row + 1
  
  writeData(wb, sheet_name, "Clinical interpretation:", startRow = current_row, startCol = 1)
  current_row <- current_row + 1
  writeData(wb, sheet_name, "  - <1.0: Most patients do not experience this event", 
            startRow = current_row, startCol = 1)
  current_row <- current_row + 1
  writeData(wb, sheet_name, "  - ~1.0: On average, each patient experiences this event once", 
            startRow = current_row, startCol = 1)
  current_row <- current_row + 1
  writeData(wb, sheet_name, "  - >1.0: Some patients experience multiple episodes of this event", 
            startRow = current_row, startCol = 1)
  current_row <- current_row + 2
  
  writeData(wb, sheet_name, ae_summary$total_incident_events, 
            startRow = current_row, startCol = 1)
  addStyle(wb, sheet_name, headerStyle, rows = current_row, 
           cols = 1:ncol(ae_summary$total_incident_events), gridExpand = TRUE)
  
  # Format numbers
  numStyle <- createStyle(numFmt = "0.0000")
  addStyle(wb, sheet_name, numStyle, 
           rows = (current_row + 1):(current_row + nrow(ae_summary$total_incident_events)), 
           cols = 2:4, gridExpand = TRUE)
  
  current_row <- current_row + nrow(ae_summary$total_incident_events) + 2
  
  # Section 2: Event Rates per 100 Patient-Years
  writeData(wb, sheet_name, "EVENT RATES PER 100 PATIENT-YEARS", 
            startRow = current_row, startCol = 1)
  addStyle(wb, sheet_name, sectionStyle, rows = current_row, cols = 1:4, gridExpand = TRUE)
  current_row <- current_row + 1
  
  writeData(wb, sheet_name, ae_summary$event_rates_per_100py, 
            startRow = current_row, startCol = 1)
  addStyle(wb, sheet_name, headerStyle, rows = current_row, 
           cols = 1:ncol(ae_summary$event_rates_per_100py), gridExpand = TRUE)
  
  # Format rates
  rateStyle <- createStyle(numFmt = "0.00")
  addStyle(wb, sheet_name, rateStyle, 
           rows = (current_row + 1):(current_row + nrow(ae_summary$event_rates_per_100py)), 
           cols = 2:4, gridExpand = TRUE)
  
  current_row <- current_row + nrow(ae_summary$event_rates_per_100py) + 1
  
  # Add patient-years information
  py_info <- paste0("Total patient-years - TCZ+GC: ", 
                    round(ae_summary$total_patient_years[1], 1),
                    " | GC alone: ", 
                    round(ae_summary$total_patient_years[2], 1))
  writeData(wb, sheet_name, py_info, startRow = current_row, startCol = 1)
  current_row <- current_row + 2
  
  # Section 3: Lifetime Per-Event Cost and QALY Impacts
  writeData(wb, sheet_name, "LIFETIME PER-EVENT COST AND QALY IMPACTS", 
            startRow = current_row, startCol = 1)
  addStyle(wb, sheet_name, sectionStyle, rows = current_row, cols = 1:7, gridExpand = TRUE)
  current_row <- current_row + 1
  
  writeData(wb, sheet_name, ae_summary$per_event_cost_qaly, 
            startRow = current_row, startCol = 1)
  addStyle(wb, sheet_name, headerStyle, rows = current_row, 
           cols = 1:ncol(ae_summary$per_event_cost_qaly), gridExpand = TRUE)
  
  # Format costs
  addStyle(wb, sheet_name, costStyle, 
           rows = (current_row + 1):(current_row + nrow(ae_summary$per_event_cost_qaly)), 
           cols = 2:4, gridExpand = TRUE)
  
  # Format QALYs (negative values)
  addStyle(wb, sheet_name, qalyStyle, 
           rows = (current_row + 1):(current_row + nrow(ae_summary$per_event_cost_qaly)), 
           cols = 5:7, gridExpand = TRUE)
  
  current_row <- current_row + nrow(ae_summary$per_event_cost_qaly) + 2
  
  # Section 4: Total Event Costs
  writeData(wb, sheet_name, "TOTAL EVENT COSTS", 
            startRow = current_row, startCol = 1)
  addStyle(wb, sheet_name, sectionStyle, rows = current_row, cols = 1:2, gridExpand = TRUE)
  current_row <- current_row + 1
  
  total_event_costs <- res$total_event_costs
  event_cost_summary <- data.frame(
    Item = c("TCZ + GC arm", "GC alone arm", "Difference", 
             "Cost savings from reduced adverse events"),
    Cost = c(
      round(total_event_costs[1], 0),
      round(total_event_costs[2], 0),
      round(total_event_costs[1] - total_event_costs[2], 0),
      round(total_event_costs[2] - total_event_costs[1], 0)
    )
  )
  
  writeData(wb, sheet_name, event_cost_summary, startRow = current_row, startCol = 1)
  addStyle(wb, sheet_name, headerStyle, rows = current_row, cols = 1:2, gridExpand = TRUE)
  
  # Format event costs
  addStyle(wb, sheet_name, createStyle(numFmt = "£#,##0"), 
           rows = (current_row + 1):(current_row + nrow(event_cost_summary)), 
           cols = 2, gridExpand = TRUE)
  
  # Set column widths
  setColWidths(wb, sheet_name, cols = 1:7, widths = "auto")
}

# Save the workbook
####################

output_file <- "03_results/scenario_analysis_adverse_events.xlsx"
saveWorkbook(wb, output_file, overwrite = TRUE)

cat("\n========================================\n")
cat("Excel file created successfully!\n")
cat(paste("Location:", output_file, "\n"))
cat("========================================\n\n")

cat("Sheets created:\n")
cat("  1. CE Summary - Cost-effectiveness summary for all scenarios\n")
cat("  2-7. AE_[scenario] - Detailed adverse event tables for each scenario\n")
cat("  8. AE Comparison - Side-by-side comparison of adverse events\n\n")

# Also save the summary table as CSV
write.csv(summary_table, "03_results/scenario_analysis_results.csv", row.names = FALSE)
cat("CSV file also saved: 03_results/scenario_analysis_results.csv\n")




#######################################
# 4. CEAC DATA SHEET
#######################################

addWorksheet(wb, "CEAC Data")

# Calculate CEAC
wtp_range <- seq(0, 50000, by = 1000)
prob_ce <- sapply(wtp_range, function(wtp) {
  mean((inc_costs_psa / inc_qalys_psa) <= wtp)
})

ceac_data <- data.frame(
  WTP_Threshold = wtp_range,
  Probability_Cost_Effective = round(prob_ce, 4)
)

writeData(wb, "CEAC Data", ceac_data, startRow = 1, startCol = 1)
addStyle(wb, "CEAC Data", headerStyle, rows = 1, cols = 1:2, gridExpand = TRUE)

#######################################
# 5. POPULATION TRACE SHEET
#######################################

addWorksheet(wb, "Population Trace")

# Create population trace for both arms
# TCZ + GC arm
pop_trace_tcz <- data.frame(
  Cycle = 0:n_cycles,
  Initial_Remission = baseline_results$pop_summary["Initial remission", , "1"],
  Relapse = baseline_results$pop_summary["Relapse", , "1"],
  Subsequent_Remission = baseline_results$pop_summary["Subsequent remission", , "1"],
  Dead = baseline_results$pop_summary["Dead", , "1"]
)

writeData(wb, "Population Trace", "TOCILIZUMAB + GC ARM", startRow = 1, startCol = 1)
writeData(wb, "Population Trace", pop_trace_tcz, startRow = 2, startCol = 1)

# GC alone arm
pop_trace_gc <- data.frame(
  Cycle = 0:n_cycles,
  Initial_Remission = baseline_results$pop_summary["Initial remission", , "0"],
  Relapse = baseline_results$pop_summary["Relapse", , "0"],
  Subsequent_Remission = baseline_results$pop_summary["Subsequent remission", , "0"],
  Dead = baseline_results$pop_summary["Dead", , "0"]
)

writeData(wb, "Population Trace", "GC ALONE ARM", startRow = n_cycles + 5, startCol = 1)
writeData(wb, "Population Trace", pop_trace_gc, startRow = n_cycles + 6, startCol = 1)

# Apply formatting
addStyle(wb, "Population Trace", headerStyle, rows = c(1, n_cycles + 5), cols = 1:5, gridExpand = TRUE)

#######################################
# 6. CYCLE-LEVEL COSTS AND QALYS
#######################################

addWorksheet(wb, "Cycle Costs and QALYs")

# TCZ + GC arm
cycle_results_tcz <- data.frame(
  Cycle = 1:n_cycles,
  State_Costs = round(baseline_results$cycle_costs[1, ], 2),
  Event_Costs = round(baseline_results$event_costs[1, ], 2),
  Total_Costs = round(baseline_results$cycle_costs[1, ] + baseline_results$event_costs[1, ], 2),
  State_QALYs = round(baseline_results$cycle_qalys[1, ], 4),
  Event_QALYs = round(baseline_results$event_qalys[1, ], 4),
  Total_QALYs = round(baseline_results$cycle_qalys[1, ] + baseline_results$event_qalys[1, ], 4)
)

writeData(wb, "Cycle Costs and QALYs", "TOCILIZUMAB + GC ARM", startRow = 1, startCol = 1)
writeData(wb, "Cycle Costs and QALYs", cycle_results_tcz, startRow = 2, startCol = 1)

# GC alone arm
cycle_results_gc <- data.frame(
  Cycle = 1:n_cycles,
  State_Costs = round(baseline_results$cycle_costs[2, ], 2),
  Event_Costs = round(baseline_results$event_costs[2, ], 2),
  Total_Costs = round(baseline_results$cycle_costs[2, ] + baseline_results$event_costs[2, ], 2),
  State_QALYs = round(baseline_results$cycle_qalys[2, ], 4),
  Event_QALYs = round(baseline_results$event_qalys[2, ], 4),
  Total_QALYs = round(baseline_results$cycle_qalys[2, ] + baseline_results$event_qalys[2, ], 4)
)

writeData(wb, "Cycle Costs and QALYs", "GC ALONE ARM", startRow = n_cycles + 5, startCol = 1)
writeData(wb, "Cycle Costs and QALYs", cycle_results_gc, startRow = n_cycles + 6, startCol = 1)

# Apply formatting
addStyle(wb, "Cycle Costs and QALYs", headerStyle, rows = c(1, n_cycles + 5), cols = 1:7, gridExpand = TRUE)

#######################################
# 7. PSA ADVERSE EVENTS DETAIL
#######################################

addWorksheet(wb, "PSA AE Details")

# Create detailed adverse event results from PSA
for (event_idx in 1:7) {
  event_name <- dimnames(psa_incident_events)$event[event_idx]
  
  psa_ae_detail <- data.frame(
    Simulation = 1:n_sim,
    TCZ_GC = round(psa_incident_events[, event_idx, 1], 4),
    GC_Alone = round(psa_incident_events[, event_idx, 2], 4),
    Difference = round(psa_incident_events[, event_idx, 1] - psa_incident_events[, event_idx, 2], 4)
  )
  
  start_row <- (event_idx - 1) * (n_sim + 3) + 1
  
  writeData(wb, "PSA AE Details", paste(event_name, "- PSA RESULTS"), 
            startRow = start_row, startCol = 1)
  writeData(wb, "PSA AE Details", psa_ae_detail, startRow = start_row + 1, startCol = 1)
  
  addStyle(wb, "PSA AE Details", headerStyle, rows = start_row, cols = 1:4, gridExpand = TRUE)
}

#######################################
# 8. COST-EFFECTIVENESS PLANE DATA
#######################################

addWorksheet(wb, "CE Plane Data")

ce_plane_data <- data.frame(
  Simulation = 1:n_sim,
  Incremental_QALYs = round(inc_qalys_psa, 4),
  Incremental_Costs = round(inc_costs_psa, 2),
  ICER = round(psa_icer, 2),
  Cost_Effective_20k = inc_costs_psa / inc_qalys_psa <= 20000,
  Cost_Effective_30k = inc_costs_psa / inc_qalys_psa <= 30000
)

writeData(wb, "CE Plane Data", ce_plane_data, startRow = 1, startCol = 1)
addStyle(wb, "CE Plane Data", headerStyle, rows = 1, cols = 1:6, gridExpand = TRUE)

#######################################
# SAVE WORKBOOK
#######################################

# Save the workbook
saveWorkbook(wb, "03_results/Model_results.xlsx", overwrite = TRUE)

cat("\n=== RESULTS EXPORTED TO EXCEL ===\n")
cat("File saved as: Model_Results.xlsx\n")
cat("Sheets created:\n")
cat("  1. Baseline Results\n")
cat("  2. PSA Summary\n")
cat("  3. PSA Detailed Results\n")
cat("  4. CEAC Data\n")
cat("  5. Population Trace\n")
cat("  6. Cycle Costs and QALYs\n")
cat("  7. PSA AE Details\n")
cat("  8. CE Plane Data\n")
