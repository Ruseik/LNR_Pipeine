################################################################################
#  02_descriptive.R — Phase 1: Descriptive Statistics & TLE Characterisation
#  Protocol Section 7.1 | SO6
#
#  Outputs: Table 1, TLE distribution, LNR distribution, outcome summary
################################################################################

source(file.path("d:/Google_SSD_RAM/_OneDrive/Research",
                 "Project_ColoRectal/LNR_20260328/Analysis/00_config.R"))

# Load analytic dataset
load(file.path(DIR_DATA, "analytic_dataset.RData"))

section_header("SCRIPT 02: DESCRIPTIVE STATISTICS (Phase 1)")

cat("Analytic cohort: N =", nrow(df), "node-positive patients\n")
cat("Analytic tier:", ANALYTIC_TIER$label, "\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 1. TABLE 1 — BASELINE CHARACTERISTICS
#    Stratified by: (a) LNR tier, (b) N-stage
#    Protocol: SMD instead of p-values for baseline balance assessment
# ═══════════════════════════════════════════════════════════════════════════════

section_header("TABLE 1: BASELINE CHARACTERISTICS")

table1_vars <- c("AG", "Sex_f", "T_stage_group", "LN_total", "LN_positive",
                 "LNR_prop", "Grade_f", "LVI_f", "TLE_adequate",
                 "Time_survival", "Outcome_OS", "Outcome_DSS")

cat_vars <- c("Sex_f", "T_stage_group", "Grade_f", "LVI_f", "TLE_adequate",
              "Outcome_OS", "Outcome_DSS")

nonnormal_vars <- c("Time_survival", "LN_total", "LN_positive", "LNR_prop")

# --- Table 1a: By LNR Tier ---
cat("--- Table 1a: Stratified by LNR Tier ---\n")
tab1_lnr <- CreateTableOne(
  vars = table1_vars,
  strata = "LNR_tier",
  data = df,
  factorVars = cat_vars,
  addOverall = TRUE
)
tab1_lnr_print <- print(tab1_lnr, 
                          nonnormal = nonnormal_vars,
                          smd = TRUE,        # SMD instead of p-values
                          test = FALSE,       # Suppress p-values per protocol
                          printToggle = TRUE,
                          showAllLevels = TRUE)

# --- Table 1b: By N-Stage ---
cat("\n--- Table 1b: Stratified by N-Stage ---\n")
tab1_nstage <- CreateTableOne(
  vars = table1_vars,
  strata = "N_stage",
  data = df,
  factorVars = cat_vars,
  addOverall = TRUE
)
tab1_nstage_print <- print(tab1_nstage,
                            nonnormal = nonnormal_vars,
                            smd = TRUE,
                            test = FALSE,
                            printToggle = TRUE,
                            showAllLevels = TRUE)

# Save tables
save_table(as.data.frame(tab1_lnr_print), "Table1_by_LNR_tier")
save_table(as.data.frame(tab1_nstage_print), "Table1_by_N_stage")

# ═══════════════════════════════════════════════════════════════════════════════
# 2. TLE CHARACTERISATION (SO6 — South Asian contextual contribution)
# ═══════════════════════════════════════════════════════════════════════════════

section_header("TLE CHARACTERISATION")

# --- 2.1 Distribution summary ---
cat("TLE Distribution (node-positive cohort):\n")
cat("  N:", sum(!is.na(df$LN_total)), "\n")
cat("  Mean (SD):", round(mean(df$LN_total, na.rm = TRUE), 1), 
    "(", round(sd(df$LN_total, na.rm = TRUE), 1), ")\n")
cat("  Median (IQR):", median(df$LN_total, na.rm = TRUE), 
    "(", quantile(df$LN_total, 0.25, na.rm = TRUE), "-",
    quantile(df$LN_total, 0.75, na.rm = TRUE), ")\n")
cat("  Range:", min(df$LN_total, na.rm = TRUE), "-", 
    max(df$LN_total, na.rm = TRUE), "\n\n")

# --- 2.2 Proportion with inadequate harvest ---
n_below_12 <- sum(df$LN_total < TLE_ADEQUATE_THRESHOLD, na.rm = TRUE)
n_below_8 <- sum(df$LN_total < TLE_STRINGENT_THRESHOLD, na.rm = TRUE)
n_tle <- sum(!is.na(df$LN_total))

ci_12 <- wilson_ci(n_below_12, n_tle)
ci_8 <- wilson_ci(n_below_8, n_tle)

cat("Harvest adequacy:\n")
cat("  TLE < 12 (AJCC minimum):  ", n_below_12, "/", n_tle, 
    " (", round(ci_12["estimate"] * 100, 1), "%; 95% CI: ",
    round(ci_12["lower"] * 100, 1), "-", round(ci_12["upper"] * 100, 1), "%)\n", sep = "")
cat("  TLE < 8 (stringent):      ", n_below_8, "/", n_tle,
    " (", round(ci_8["estimate"] * 100, 1), "%; 95% CI: ",
    round(ci_8["lower"] * 100, 1), "-", round(ci_8["upper"] * 100, 1), "%)\n", sep = "")

# --- 2.3 TLE histogram ---
p_tle_hist <- ggplot(df, aes(x = LN_total)) +
  geom_histogram(binwidth = 2, fill = "#4393C3", colour = "white", alpha = 0.85) +
  geom_vline(xintercept = TLE_ADEQUATE_THRESHOLD, linetype = "dashed", 
             colour = "#D95F02", linewidth = 0.8) +
  annotate("text", x = TLE_ADEQUATE_THRESHOLD + 1, 
           y = Inf, vjust = 2, hjust = 0,
           label = paste0("AJCC threshold = ", TLE_ADEQUATE_THRESHOLD),
           colour = "#D95F02", fontface = "italic", size = 3.5) +
  labs(title = "Distribution of Total Lymph Nodes Examined (TLE)",
       subtitle = "Node-positive CRC cohort",
       x = "Total Lymph Nodes Examined",
       y = "Number of Patients") +
  theme_publication()
save_figure(p_tle_hist, "Fig_TLE_histogram")

# --- 2.4 TLE by N-stage ---
p_tle_nstage <- ggplot(df, aes(x = N_stage, y = LN_total, fill = N_stage)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  geom_hline(yintercept = TLE_ADEQUATE_THRESHOLD, linetype = "dashed",
             colour = "#D95F02", linewidth = 0.5) +
  scale_fill_manual(values = COL_NSTAGE) +
  labs(title = "Total Lymph Nodes Examined by N-Stage",
       x = "N-Stage", y = "Total Lymph Nodes Examined") +
  guides(fill = "none") +
  theme_publication()
save_figure(p_tle_nstage, "Fig_TLE_by_Nstage")

# --- 2.5 TLE by LNR tier ---
p_tle_lnr <- ggplot(df, aes(x = LNR_tier, y = LN_total, fill = LNR_tier)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  geom_hline(yintercept = TLE_ADEQUATE_THRESHOLD, linetype = "dashed",
             colour = "#D95F02", linewidth = 0.5) +
  scale_fill_manual(values = COL_LNR_TIERS) +
  labs(title = "Total Lymph Nodes Examined by LNR Tier",
       x = "LNR Tier (Rosenberg)", y = "Total Lymph Nodes Examined") +
  guides(fill = "none") +
  theme_publication()
save_figure(p_tle_lnr, "Fig_TLE_by_LNR_tier")

# ═══════════════════════════════════════════════════════════════════════════════
# 3. LNR DISTRIBUTION
# ═══════════════════════════════════════════════════════════════════════════════

section_header("LNR DISTRIBUTION")

cat("LNR (proportion) summary:\n")
cat("  Mean (SD):", round(mean(df$LNR_prop, na.rm = TRUE), 4), 
    "(", round(sd(df$LNR_prop, na.rm = TRUE), 4), ")\n")
cat("  Median (IQR):", round(median(df$LNR_prop, na.rm = TRUE), 4),
    "(", round(quantile(df$LNR_prop, 0.25, na.rm = TRUE), 4), "-",
    round(quantile(df$LNR_prop, 0.75, na.rm = TRUE), 4), ")\n")
cat("  Range:", round(min(df$LNR_prop, na.rm = TRUE), 4), "-",
    round(max(df$LNR_prop, na.rm = TRUE), 4), "\n\n")

cat("LNR Tier Distribution (Rosenberg 2008):\n")
tier_tab <- as.data.frame(table(df$LNR_tier))
names(tier_tab) <- c("LNR_Tier", "N")
tier_tab$Pct <- round(100 * tier_tab$N / sum(tier_tab$N), 1)
print(tier_tab)
save_table(tier_tab, "LNR_tier_distribution")

# --- LNR histogram ---
p_lnr_hist <- ggplot(df, aes(x = LNR_prop)) +
  geom_histogram(binwidth = 0.05, fill = "#B2182B", colour = "white", alpha = 0.85) +
  geom_vline(xintercept = LNR_CUT_LOW_INTER, linetype = "dashed", 
             colour = "#4DAF4A", linewidth = 0.8) +
  geom_vline(xintercept = LNR_CUT_INTER_HIGH, linetype = "dashed",
             colour = "#FF7F00", linewidth = 0.8) +
  annotate("text", x = LNR_CUT_LOW_INTER + 0.01, y = Inf, vjust = 2, hjust = 0,
           label = "0.05", colour = "#4DAF4A", fontface = "italic", size = 3.5) +
  annotate("text", x = LNR_CUT_INTER_HIGH + 0.01, y = Inf, vjust = 2, hjust = 0,
           label = "0.20", colour = "#FF7F00", fontface = "italic", size = 3.5) +
  labs(title = "Distribution of Lymph Node Ratio (LNR)",
       subtitle = "Node-positive CRC cohort | Rosenberg thresholds shown",
       x = "LNR (proportion)",
       y = "Number of Patients") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme_publication()
save_figure(p_lnr_hist, "Fig_LNR_histogram")

# ═══════════════════════════════════════════════════════════════════════════════
# 4. OUTCOME SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════

section_header("OUTCOME SUMMARY")

n_total <- nrow(df)
n_os <- sum(df$Outcome_OS == 1, na.rm = TRUE)
n_dss <- sum(df$Outcome_DSS == 1, na.rm = TRUE)

cat("Overall Survival:\n")
cat("  Total patients:", n_total, "\n")
cat("  Deaths (any cause):", n_os, "(", round(100 * n_os / n_total, 1), "%)\n")
cat("  Censored:", n_total - n_os, "(", round(100 * (n_total - n_os) / n_total, 1), "%)\n\n")

cat("Disease-Specific Survival:\n")
cat("  Cancer deaths:", n_dss, "(", round(100 * n_dss / n_total, 1), "%)\n")
cat("  Non-cancer deaths:", n_os - n_dss, "\n")

# -- CSS data completeness check --
if (n_os > 0) {
  css_completeness <- n_dss / n_os  # Proportion of deaths with cancer-specific attribution
  # Note: this is an approximation — assumes DSS == 1 implies cancer death
  cat("\n  CSS data completeness proxy:", round(css_completeness * 100, 1), "%\n")
  cat("  (Proportion of deceased patients with DSS event recorded)\n")
  if (css_completeness >= CSS_COD_THRESHOLD) {
    cat("  --> CSS analysis PERMITTED (>= 85% threshold)\n")
  } else {
    cat("  --> CSS analysis may be limited (< 85% threshold)\n")
  }
}

# -- Median follow-up (reverse KM method) --
surv_fu <- survfit(Surv(Time_survival, 1 - Outcome_OS) ~ 1, data = df)
median_fu <- summary(surv_fu)$table["median"]
cat("\nMedian follow-up (reverse KM method):", 
    round(median_fu, 0), "days (",
    round(median_fu / 365.25, 1), "years)\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 5. CROSS-TABULATION: N-STAGE × LNR TIER (Preview for reclassification)
# ═══════════════════════════════════════════════════════════════════════════════

section_header("N-STAGE × LNR TIER CROSS-TABULATION")

cross_tab <- table(df$N_stage, df$LNR_tier)
cat("Counts:\n")
print(cross_tab)
cat("\nRow percentages (within each N-stage):\n")
print(round(prop.table(cross_tab, margin = 1) * 100, 1))

save_table(as.data.frame.matrix(cross_tab), "Nstage_LNRtier_crosstab")

# ═══════════════════════════════════════════════════════════════════════════════
# 6. PATIENT FLOW DIAGRAM (Formatted Text)
# ═══════════════════════════════════════════════════════════════════════════════

section_header("CONSORT-STYLE PATIENT FLOW")

cat("┌─────────────────────────────────────────────────┐\n")
cat("│  Total records in database: ", sprintf("%-20d", flow$step1_total), "│\n")
cat("└──────────────────────┬──────────────────────────┘\n")
cat("                       │\n")
cat("  Excluded: invalid LN │ n =", flow$excluded_invalid_ln, "\n")
cat("                       ▼\n")
cat("┌─────────────────────────────────────────────────┐\n")
cat("│  Valid LN records:    ", sprintf("%-20d", flow$step2_after_ln_valid), "│\n")
cat("└──────────────────────┬──────────────────────────┘\n")
cat("                       │\n")
cat("  Excluded: missing    │ n =", flow$excluded_missing_surv, "\n")
cat("  survival data        │\n")
cat("                       ▼\n")
cat("┌─────────────────────────────────────────────────┐\n")
cat("│  Complete survival:   ", sprintf("%-20d", flow$step3_after_surv), "│\n")
cat("└──────────────────────┬──────────────────────────┘\n")
cat("                       │\n")
cat("  Excluded: node-neg   │ n =", flow$excluded_node_negative, "\n")
cat("  (", NODE_POSITIVE_FILTER, ")             │\n", sep = "")
cat("                       ▼\n")
cat("┌─────────────────────────────────────────────────┐\n")
cat("│  Node-positive:      ", sprintf("%-20d", flow$step4_node_positive), "│\n")
cat("└──────────────────────┬──────────────────────────┘\n")
cat("                       │\n")
cat("  Excluded: TLE missing│ n =", flow$excluded_tle_missing, "\n")
cat("                       ▼\n")
cat("┌═════════════════════════════════════════════════┐\n")
cat("│  FINAL ANALYTIC COHORT: ", sprintf("%-18d", flow$step5_final_analytic), "│\n")
cat("│  OS events: ", sprintf("%-8d", n_os), 
    "  DSS events: ", sprintf("%-8d", n_dss), "│\n")
cat("└═════════════════════════════════════════════════┘\n")

cat("\n", strrep("=", 78), "\n")
cat("  SCRIPT 02 COMPLETE\n")
cat(strrep("=", 78), "\n")
