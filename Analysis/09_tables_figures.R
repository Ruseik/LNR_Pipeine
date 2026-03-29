################################################################################
#  09_tables_figures.R — Publication Tables & Figures Compilation
#  Assembles all outputs into final publication-ready format.
#
#  Run AFTER all preceding scripts (01–08) have completed successfully.
################################################################################

source(file.path("d:/Google_SSD_RAM/_OneDrive/Research",
                 "Project_ColoRectal/LNR_20260328/Analysis/00_config.R"))

load(file.path(DIR_DATA, "analytic_dataset.RData"))

# Load all model results (wrapped in tryCatch for robustness)
tryCatch(load(file.path(DIR_MODELS, "phase2_unadjusted_results.RData")),
         error = function(e) cat("Phase 2 results not found.\n"))
tryCatch(load(file.path(DIR_MODELS, "phase3_cox_models.RData")),
         error = function(e) cat("Phase 3 model results not found.\n"))
tryCatch(load(file.path(DIR_MODELS, "phase3_performance_metrics.RData")),
         error = function(e) cat("Phase 3 performance results not found.\n"))
tryCatch(load(file.path(DIR_MODELS, "phase4_reclassification.RData")),
         error = function(e) cat("Phase 4 results not found.\n"))
tryCatch(load(file.path(DIR_MODELS, "phase5_subgroup_sensitivity.RData")),
         error = function(e) cat("Phase 5 results not found.\n"))
tryCatch(load(file.path(DIR_MODELS, "phase_css_competing_risks.RData")),
         error = function(e) cat("CSS results not found.\n"))

section_header("SCRIPT 09: TABLES & FIGURES COMPILATION")

# ═══════════════════════════════════════════════════════════════════════════════
# 1. MANUSCRIPT FIGURE 1: CONSORT FLOW DIAGRAM (Text-based)
# ═══════════════════════════════════════════════════════════════════════════════

section_header("FIGURE 1: CONSORT FLOW")

flow_text <- paste0(
  "CONSORT-Style Patient Flow Diagram\n",
  "===================================\n\n",
  "Total records in database: ", flow$step1_total, "\n",
  "  └─ Excluded (invalid LN data): ", flow$excluded_invalid_ln, "\n",
  "  └─ Excluded (missing survival): ", flow$excluded_missing_surv, "\n",
  "  └─ Excluded (node-negative): ", flow$excluded_node_negative, "\n",
  "     (Filter: ", NODE_POSITIVE_FILTER, ")\n",
  "  └─ Excluded (TLE missing/zero): ", flow$excluded_tle_missing, "\n",
  "  ────────────────────\n",
  "  FINAL ANALYTIC COHORT: ", flow$step5_final_analytic, "\n",
  "  OS events: ", sum(df$Outcome_OS == 1), "\n",
  "  DSS events: ", sum(df$Outcome_DSS == 1), "\n",
  "  Analytic tier: ", ANALYTIC_TIER$label, "\n"
)
cat(flow_text)
writeLines(flow_text, file.path(DIR_TABLES, "Figure1_CONSORT_flow.txt"))

# ═══════════════════════════════════════════════════════════════════════════════
# 2. MANUSCRIPT TABLE 2: MULTIVARIABLE COX MODEL RESULTS
# ═══════════════════════════════════════════════════════════════════════════════

section_header("TABLE 2: MULTIVARIABLE COX MODEL RESULTS")

if (exists("models") && !is.null(models$C) && !is.null(models$D)) {
  # Model C (N-stage full MVA)
  cox_C_tidy <- coxph(Surv(Time_survival, Outcome_OS) ~ 
                        N_stage + LN_total + AG + Sex_f + 
                        T_stage_group + Grade_f + LVI_f, data = df)
  tidy_C <- broom::tidy(cox_C_tidy, exponentiate = TRUE, conf.int = TRUE)
  tidy_C$Model <- "C (N-stage)"
  
  # Model D (LNR full MVA)
  cox_D_tidy <- coxph(Surv(Time_survival, Outcome_OS) ~ 
                        LNR_prop + LN_total + AG + Sex_f + 
                        T_stage_group + Grade_f + LVI_f, data = df)
  tidy_D <- broom::tidy(cox_D_tidy, exponentiate = TRUE, conf.int = TRUE)
  tidy_D$Model <- "D (LNR)"
  
  table2 <- rbind(tidy_C, tidy_D)
  table2$HR_CI <- sprintf("%.2f (%.2f–%.2f)", 
                           table2$estimate, table2$conf.low, table2$conf.high)
  table2$p_formatted <- sapply(table2$p.value, format_pval)
  
  table2_export <- table2[, c("Model", "term", "HR_CI", "p_formatted")]
  names(table2_export) <- c("Model", "Variable", "HR (95% CI)", "p-value")
  
  cat("Table 2: Full MVA Cox Model Results\n")
  print(table2_export)
  save_table(table2_export, "Table2_MVA_results")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 3. MANUSCRIPT TABLE 3: C-STATISTIC COMPARISON ACROSS LAYERS
# ═══════════════════════════════════════════════════════════════════════════════

section_header("TABLE 3: C-STATISTIC COMPARISON")

if (exists("c_results") && nrow(c_results) > 0) {
  table3 <- c_results[, c("Model", "C_statistic", "CI_lower", "CI_upper")]
  table3$C_CI <- sprintf("%.4f (%.4f–%.4f)", 
                          table3$C_statistic, table3$CI_lower, table3$CI_upper)
  
  # Add layer info
  table3$Layer <- ifelse(table3$Model %in% c("A", "B"), "1 (Unadjusted)",
                  ifelse(table3$Model %in% c("C2", "D2"), "2 (TLE-adjusted)",
                         "3 (Full MVA)"))
  
  # Add ΔC within each layer
  table3$Delta_C <- NA
  pairs <- list(c("A", "B"), c("C2", "D2"), c("C", "D"))
  for (pair in pairs) {
    if (all(pair %in% table3$Model)) {
      c1 <- table3$C_statistic[table3$Model == pair[1]]
      c2 <- table3$C_statistic[table3$Model == pair[2]]
      table3$Delta_C[table3$Model == pair[2]] <- round(c2 - c1, 4)
    }
  }
  
  table3_export <- table3[, c("Layer", "Model", "C_CI", "Delta_C")]
  names(table3_export) <- c("Layer", "Model", "C-statistic (95% CI)", 
                             "ΔC (LNR - N-stage)")
  
  cat("Table 3: C-Statistic Comparison\n")
  print(table3_export)
  save_table(table3_export, "Table3_C_statistic_comparison")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 4. MANUSCRIPT TABLE 4: RECLASSIFICATION TABLE
# ═══════════════════════════════════════════════════════════════════════════════

section_header("TABLE 4: RECLASSIFICATION")

if (exists("cross_tab")) {
  cat("N-Stage × LNR Tier Cross-Tabulation:\n")
  print(addmargins(cross_tab))
  
  reclass_output <- as.data.frame.matrix(addmargins(cross_tab))
  save_table(reclass_output, "Table4_reclassification")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 5. MULTI-PANEL COMPOSITE FIGURES
# ═══════════════════════════════════════════════════════════════════════════════

section_header("COMPOSITE FIGURES")

# Try to create a composite KM figure (N-stage + LNR tier side by side)
tryCatch({
  km_ns <- survfit(Surv(Time_survival, Outcome_OS) ~ N_stage, data = df)
  km_lnr <- survfit(Surv(Time_survival, Outcome_OS) ~ LNR_tier, data = df)
  
  p1 <- ggsurvplot(km_ns, data = df, palette = unname(COL_NSTAGE),
                    pval = TRUE, conf.int = TRUE, conf.int.alpha = 0.1,
                    xlab = "Time (days)", ylab = "OS Probability",
                    title = "A. Overall Survival by N-Stage",
                    legend.title = "N-Stage", legend.labs = c("N1", "N2"),
                    ggtheme = theme_publication(base_size = 10),
                    break.time.by = 365.25, legend = "bottom")
  
  p2 <- ggsurvplot(km_lnr, data = df, palette = unname(COL_LNR_TIERS),
                    pval = TRUE, conf.int = TRUE, conf.int.alpha = 0.1,
                    xlab = "Time (days)", ylab = "OS Probability",
                    title = "B. Overall Survival by LNR Tier",
                    legend.title = "LNR Tier",
                    ggtheme = theme_publication(base_size = 10),
                    break.time.by = 365.25, legend = "bottom")
  
  # Combine using cowplot or direct arrangement
  pdf(file.path(DIR_FIGURES, "Fig_composite_KM_OS.pdf"),
      width = 14, height = 7)
  par(mfrow = c(1, 2))
  # Since ggsurvplot returns a list, we print them separately
  gridExtra::grid.arrange(p1$plot, p2$plot, ncol = 2)
  dev.off()
  
  cat("Saved: Fig_composite_KM_OS (combined N-stage + LNR tier)\n")
}, error = function(e) {
  cat("Composite figure generation failed:", e$message, "\n")
})

# ═══════════════════════════════════════════════════════════════════════════════
# 6. COMPREHENSIVE RESULTS SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════

section_header("COMPREHENSIVE RESULTS SUMMARY")

cat("═══════════════════════════════════════════════════════════\n")
cat("  LNR vs N-Stage Prognostic Comparison — Key Findings\n")
cat("═══════════════════════════════════════════════════════════\n\n")

cat("Study Population:\n")
cat("  Analytic cohort: N =", nrow(df), "node-positive CRC patients\n")
cat("  Node-positive filter:", NODE_POSITIVE_FILTER, "\n")
cat("  OS events:", sum(df$Outcome_OS == 1), "\n")
cat("  DSS events:", sum(df$Outcome_DSS == 1), "\n")
cat("  Analytic tier:", ANALYTIC_TIER$label, "\n\n")

cat("TLE Characterisation (South Asian contextual data):\n")
cat("  Median TLE:", median(df$LN_total, na.rm = TRUE), "nodes\n")
cat("  TLE < 12 (AJCC minimum):", 
    round(100 * sum(df$LN_total < 12) / nrow(df), 1), "%\n")
cat("  TLE < 8 (stringent):", 
    round(100 * sum(df$LN_total < 8) / nrow(df), 1), "%\n\n")

cat("LNR Distribution:\n")
cat("  Median LNR:", round(median(df$LNR_prop, na.rm = TRUE), 3), "\n")
cat("  LNR tier distribution:\n")
print(table(df$LNR_tier))

if (exists("c_results") && nrow(c_results) > 0) {
  cat("\nC-Statistic Summary:\n")
  for (i in seq_len(nrow(c_results))) {
    cat(sprintf("  Model %s: C = %.4f\n", c_results$Model[i], 
                c_results$C_statistic[i]))
  }
}

cat("\n═══════════════════════════════════════════════════════════\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 7. FILE INVENTORY
# ═══════════════════════════════════════════════════════════════════════════════

section_header("OUTPUT FILE INVENTORY")

cat("Tables:\n")
table_files <- list.files(DIR_TABLES, full.names = FALSE)
for (f in table_files) cat("  ", f, "\n")

cat("\nFigures:\n")
fig_files <- list.files(DIR_FIGURES, full.names = FALSE)
for (f in fig_files) cat("  ", f, "\n")

cat("\nModels:\n")
model_files <- list.files(DIR_MODELS, full.names = FALSE)
for (f in model_files) cat("  ", f, "\n")

cat("\nAppendix:\n")
app_files <- list.files(DIR_APPENDIX, full.names = FALSE)
for (f in app_files) cat("  ", f, "\n")

cat("\n", strrep("=", 78), "\n")
cat("  SCRIPT 09 COMPLETE — ALL OUTPUTS COMPILED\n")
cat("  Total tables:", length(table_files), "\n")
cat("  Total figures:", length(fig_files), "\n")
cat(strrep("=", 78), "\n")
