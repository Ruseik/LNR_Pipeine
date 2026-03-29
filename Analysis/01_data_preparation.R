################################################################################
#  01_data_preparation.R — Data Loading, Cleaning, and Cohort Selection
#  Protocol: LNR vs TNM N-Stage Prognostic Comparison in Node-Positive CRC
#
#  Inputs:  Raw CSV (CRCALVS_8_1.csv)
#  Outputs: Cleaned analytic dataset (.RData + .csv), CONSORT flow, gating tier
#
#  Source 00_config.R before running this script.
################################################################################

source(file.path("d:/Google_SSD_RAM/_OneDrive/Research",
                 "Project_ColoRectal/LNR_20260328/Analysis/00_config.R"))

section_header("SCRIPT 01: DATA PREPARATION")

# ═══════════════════════════════════════════════════════════════════════════════
# 1. LOAD RAW DATA
# ═══════════════════════════════════════════════════════════════════════════════

cat("Loading raw data from:", FILE_RAW_DATA, "\n")
df_raw <- read.csv(FILE_RAW_DATA, header = TRUE, stringsAsFactors = FALSE)

cat("  Raw dataset dimensions:", nrow(df_raw), "rows x", ncol(df_raw), "columns\n")

# Drop useless column X
if ("X" %in% names(df_raw)) {
  df_raw$X <- NULL
  cat("  Dropped column 'X' (empty/useless).\n")
}

# Remove any completely empty trailing rows
df_raw <- df_raw[!is.na(df_raw$Index) & df_raw$Index != "", ]
cat("  After removing empty rows:", nrow(df_raw), "patients\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 2. DATA INTEGRITY CHECKS
# ═══════════════════════════════════════════════════════════════════════════════

section_header("DATA INTEGRITY CHECKS")

# --- 2.1 Check for duplicate Index values ---
n_dup <- sum(duplicated(df_raw$Index))
cat("Duplicate patient IDs:", n_dup, "\n")
if (n_dup > 0) {
  cat("  WARNING: Duplicates found at Index:", 
      df_raw$Index[duplicated(df_raw$Index)], "\n")
}

# --- 2.2 Check LN_positive <= LN_total ---
invalid_ln <- df_raw$LN_positive > df_raw$LN_total
n_invalid_ln <- sum(invalid_ln, na.rm = TRUE)
cat("Records where LN_positive > LN_total:", n_invalid_ln, "\n")
if (n_invalid_ln > 0) {
  cat("  *** DATA ERROR: The following records have positive > total nodes:\n")
  print(df_raw[which(invalid_ln), c("Index", "LN_total", "LN_positive", "LNR")])
  cat("  These records will be EXCLUDED per protocol (Section 4.4).\n")
}

# --- 2.3 Validate LNR against recomputed value ---
# Dataset LNR is stored as percentage (0-100)
df_raw$LNR_recomputed <- ifelse(df_raw$LN_total > 0,
                                 (df_raw$LN_positive / df_raw$LN_total) * 100,
                                 0)
lnr_discrepancy <- abs(df_raw$LNR - df_raw$LNR_recomputed) > 0.5
n_discrepancy <- sum(lnr_discrepancy, na.rm = TRUE)
cat("LNR discrepancies (|stored - recomputed| > 0.5%):", n_discrepancy, "\n")
if (n_discrepancy > 0) {
  cat("  Records with discrepancies:\n")
  print(df_raw[which(lnr_discrepancy), 
               c("Index", "LN_total", "LN_positive", "LNR", "LNR_recomputed")])
}

# --- 2.4 Check value ranges ---
cat("\nValue range checks:\n")
cat("  Time_survival: [", min(df_raw$Time_survival, na.rm = TRUE), ",", 
    max(df_raw$Time_survival, na.rm = TRUE), "] days\n")
cat("  Age (AG): [", min(df_raw$AG, na.rm = TRUE), ",", 
    max(df_raw$AG, na.rm = TRUE), "] years\n")
cat("  T-stage (ST_main): ", paste(sort(unique(df_raw$ST_main)), collapse = ", "), "\n")
cat("  N-stage (SN_main): ", paste(sort(unique(df_raw$SN_main)), collapse = ", "), "\n")
cat("  LN_total: [", min(df_raw$LN_total, na.rm = TRUE), ",", 
    max(df_raw$LN_total, na.rm = TRUE), "]\n")
cat("  LN_positive: [", min(df_raw$LN_positive, na.rm = TRUE), ",", 
    max(df_raw$LN_positive, na.rm = TRUE), "]\n")
cat("  Grade: ", paste(sort(unique(df_raw$Grade)), collapse = ", "), "\n")
cat("  LVI: ", paste(sort(unique(df_raw$LVI)), collapse = ", "), "\n")
cat("  Outcome_OS: ", paste(sort(unique(df_raw$Outcome_OS)), collapse = ", "), "\n")
cat("  Outcome_DSS: ", paste(sort(unique(df_raw$Outcome_DSS)), collapse = ", "), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 3. MISSINGNESS ASSESSMENT
# ═══════════════════════════════════════════════════════════════════════════════

section_header("MISSINGNESS ASSESSMENT")

analysis_vars <- c("Time_survival", "Outcome_OS", "Outcome_DSS", 
                    "AG", "GN", "ST_main", "SN_main", 
                    "LN_total", "LN_positive", "LNR", "Grade", "LVI")

missing_summary <- data.frame(
  Variable = analysis_vars,
  N_missing = sapply(analysis_vars, function(v) sum(is.na(df_raw[[v]]))),
  Pct_missing = sapply(analysis_vars, function(v) 
    round(100 * sum(is.na(df_raw[[v]])) / nrow(df_raw), 2))
)
print(missing_summary)

any_high_missing <- any(missing_summary$Pct_missing > MICE_MISSINGNESS_THRESHOLD * 100)
cat("\nAny variable with >", MICE_MISSINGNESS_THRESHOLD * 100, "% missingness:", 
    any_high_missing, "\n")
if (any_high_missing) {
  cat("  --> MICE multiple imputation will be applied in downstream scripts.\n")
} else {
  cat("  --> Complete case analysis is the primary approach (no imputation needed).\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 4. CONSORT-STYLE PATIENT FLOW (Tracking all exclusions)
# ═══════════════════════════════════════════════════════════════════════════════

section_header("CONSORT-STYLE PATIENT FLOW")

flow <- list()
flow$step1_total <- nrow(df_raw)
cat("Step 1 — Total records in database:", flow$step1_total, "\n")

# --- 4.1 Exclude records with invalid LN data ---
df_work <- df_raw[!invalid_ln | is.na(invalid_ln), ]
flow$excluded_invalid_ln <- flow$step1_total - nrow(df_work)
flow$step2_after_ln_valid <- nrow(df_work)
cat("Step 2 — After excluding invalid LN (positive > total):", 
    flow$step2_after_ln_valid, 
    "(excluded:", flow$excluded_invalid_ln, ")\n")

# --- 4.2 Exclude records with missing survival data ---
has_surv <- !is.na(df_work$Time_survival) & !is.na(df_work$Outcome_OS)
df_work <- df_work[has_surv, ]
flow$excluded_missing_surv <- flow$step2_after_ln_valid - nrow(df_work)
flow$step3_after_surv <- nrow(df_work)
cat("Step 3 — After excluding missing survival data:", 
    flow$step3_after_surv,
    "(excluded:", flow$excluded_missing_surv, ")\n")

# --- 4.3 Store the FULL dataset (node-positive + node-negative) ---
df_all <- df_work
cat("\n  Full cohort (all nodes) retained for N-stage-only models: N =", 
    nrow(df_all), "\n")

# --- 4.4 Filter to NODE-POSITIVE cohort ---
cat("\n  >>> Filtering strategy: ", NODE_POSITIVE_FILTER, " <<<\n")

if (NODE_POSITIVE_FILTER == "SN_main") {
  node_pos_mask <- df_work$SN_main >= 1
  filter_description <- "SN_main >= 1 (N-stage N1 or N2)"
} else if (NODE_POSITIVE_FILTER == "LN_positive") {
  node_pos_mask <- df_work$LN_positive >= 1
  filter_description <- "LN_positive >= 1 (at least 1 positive lymph node)"
} else {
  stop("Invalid NODE_POSITIVE_FILTER value. Must be 'SN_main' or 'LN_positive'.")
}

df_nodepos <- df_work[node_pos_mask, ]
flow$excluded_node_negative <- nrow(df_work) - nrow(df_nodepos)
flow$step4_node_positive <- nrow(df_nodepos)
cat("Step 4 — Node-positive cohort (", filter_description, "):", 
    flow$step4_node_positive,
    "(excluded node-negative:", flow$excluded_node_negative, ")\n")

# --- 4.5 Exclude patients with missing TLE (cannot compute LNR) ---
tle_missing_mask <- is.na(df_nodepos$LN_total) | df_nodepos$LN_total == 0
n_tle_missing <- sum(tle_missing_mask)
if (n_tle_missing > 0) {
  cat("  NOTE:", n_tle_missing, "node-positive patients have missing/zero TLE. ",
      "They are excluded from LNR models but retained for N-stage-only.\n")
  df_nodepos_tle_avail <- df_nodepos[!tle_missing_mask, ]
} else {
  df_nodepos_tle_avail <- df_nodepos
}
flow$excluded_tle_missing <- n_tle_missing
flow$step5_final_analytic <- nrow(df_nodepos_tle_avail)
cat("Step 5 — Final analytic cohort (LNR computable):", 
    flow$step5_final_analytic, "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 5. CREATE DERIVED VARIABLES
# ═══════════════════════════════════════════════════════════════════════════════

section_header("CREATING DERIVED VARIABLES")

df <- df_nodepos_tle_avail  # Working analytic dataset

# --- 5.1 LNR as proportion (0-1) — recomputed from raw counts ---
df$LNR_prop <- df$LN_positive / df$LN_total
cat("LNR_prop (proportion 0-1) created. Range: [", 
    round(min(df$LNR_prop, na.rm = TRUE), 4), ",", 
    round(max(df$LNR_prop, na.rm = TRUE), 4), "]\n")

# --- 5.2 LNR categorical tiers (Rosenberg 2008) ---
df$LNR_tier <- cut(df$LNR_prop, 
                    breaks = c(-Inf, LNR_CUT_LOW_INTER, LNR_CUT_INTER_HIGH, Inf),
                    labels = LNR_TIER_LABELS,
                    right = FALSE)
cat("LNR tier distribution:\n")
print(table(df$LNR_tier, useNA = "ifany"))

# --- 5.3 Alternative LNR categorisations for sensitivity ---
df$LNR_tier_alt1 <- cut(df$LNR_prop,
                         breaks = c(-Inf, LNR_ALT_CUT_1, Inf),
                         labels = c("Low (<0.05)", "Intermediate (0.05-0.40)", 
                                    "High (>0.40)"),
                         right = FALSE)

df$LNR_tier_binary <- ifelse(df$LNR_prop < LNR_ALT_CUT_2, 
                              "Low (<0.20)", "High (>=0.20)")
df$LNR_tier_binary <- factor(df$LNR_tier_binary, 
                              levels = c("Low (<0.20)", "High (>=0.20)"))

# --- 5.4 N-stage factor ---
df$N_stage <- factor(df$SN_main, levels = c(1, 2), labels = c("N1", "N2"))
cat("\nN-stage distribution:\n")
print(table(df$N_stage, useNA = "ifany"))

# --- 5.5 T-stage grouped (per protocol: T1/T2 vs T3 vs T4) ---
df$T_stage_group <- factor(
  ifelse(df$ST_main <= 2, "T1/T2",
         ifelse(df$ST_main == 3, "T3", "T4")),
  levels = c("T1/T2", "T3", "T4")
)
cat("\nT-stage group distribution:\n")
print(table(df$T_stage_group, useNA = "ifany"))

# --- 5.6 TLE adequacy ---
df$TLE_adequate <- ifelse(df$LN_total >= TLE_ADEQUATE_THRESHOLD, 
                           "Adequate (>=12)", "Inadequate (<12)")
df$TLE_adequate <- factor(df$TLE_adequate, 
                           levels = c("Inadequate (<12)", "Adequate (>=12)"))
df$TLE_adequate_binary <- as.integer(df$LN_total >= TLE_ADEQUATE_THRESHOLD)
cat("\nTLE adequacy distribution:\n")
print(table(df$TLE_adequate, useNA = "ifany"))

# Stringent threshold
df$TLE_stringent <- ifelse(df$LN_total >= TLE_STRINGENT_THRESHOLD,
                            "Adequate (>=8)", "Inadequate (<8)")
df$TLE_stringent <- factor(df$TLE_stringent,
                            levels = c("Inadequate (<8)", "Adequate (>=8)"))

# --- 5.7 Grade factor ---
df$Grade_f <- factor(df$Grade, levels = c(0, 1), 
                      labels = c("Low grade", "High grade"))
if (all(df$Grade %in% c(1, 2), na.rm = TRUE)) {
  # Re-check: data shows Grade = 1 (low/high? 2=high grade)
  # Based on user: 0=low grade, 1=high grade
  # But data shows 1 and 2 — need to re-map
  df$Grade_f <- factor(df$Grade, levels = c(1, 2), 
                        labels = c("Low grade", "High grade"))
}
cat("\nGrade distribution:\n")
print(table(df$Grade_f, useNA = "ifany"))

# --- 5.8 LVI factor ---
df$LVI_f <- factor(df$LVI, levels = c(0, 1), 
                    labels = c("Absent", "Present"))
cat("\nLVI distribution:\n")
print(table(df$LVI_f, useNA = "ifany"))

# --- 5.9 Sex factor ---
df$Sex_f <- factor(df$GN, levels = c(0, 1), labels = c("Male", "Female"))
cat("\nSex distribution:\n")
print(table(df$Sex_f, useNA = "ifany"))

# --- 5.10 Survival time in months (for display) ---
df$Time_survival_months <- df$Time_survival / 30.44
df$Time_survival_years <- df$Time_survival / 365.25

# ═══════════════════════════════════════════════════════════════════════════════
# 6. EVENT-COUNT GATING
# ═══════════════════════════════════════════════════════════════════════════════

section_header("EVENT-COUNT GATING")

n_os_events <- sum(df$Outcome_OS == 1, na.rm = TRUE)
n_dss_events <- sum(df$Outcome_DSS == 1, na.rm = TRUE)

cat("Node-positive analytic cohort: N =", nrow(df), "\n")
cat("Overall Survival events (deaths):", n_os_events, "\n")
cat("Disease-Specific Survival events:", n_dss_events, "\n")
cat("Median follow-up (days):", round(median(df$Time_survival), 0), "\n")
cat("Median follow-up (years):", round(median(df$Time_survival_years), 2), "\n\n")

ANALYTIC_TIER <- get_analytic_tier(n_os_events)
cat("═══ ASSIGNED ANALYTIC TIER ═══\n")
cat("  Tier:", ANALYTIC_TIER$label, "\n")
cat("  Events:", ANALYTIC_TIER$n_events, "\n")
cat("  Max degrees of freedom:", ANALYTIC_TIER$max_df, "\n")
cat("  Permitted: MVA =", ANALYTIC_TIER$permit_mva, 
    "| RCS =", ANALYTIC_TIER$permit_rcs,
    "| DCA =", ANALYTIC_TIER$permit_dca,
    "| Interaction =", ANALYTIC_TIER$permit_interaction, "\n")
cat("  Description:", ANALYTIC_TIER$description, "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 7. SUMMARY OF ANALYTIC DATASET
# ═══════════════════════════════════════════════════════════════════════════════

section_header("ANALYTIC DATASET SUMMARY")

cat("Final analytic dataset: N =", nrow(df), "node-positive patients\n\n")

cat("Variable listing:\n")
str(df[, c("Index", "Time_survival", "Outcome_OS", "Outcome_DSS",
           "AG", "Sex_f", "T_stage_group", "N_stage", 
           "LN_total", "LN_positive", "LNR_prop", "LNR_tier",
           "Grade_f", "LVI_f", "TLE_adequate")])

# ═══════════════════════════════════════════════════════════════════════════════
# 8. IDENTIFY DISCORDANCES BETWEEN FILTER METHODS
# ═══════════════════════════════════════════════════════════════════════════════

section_header("FILTER METHOD COMPARISON")

# Show patients where SN_main and LN_positive disagree on node-positivity
discordant <- df_work[(df_work$SN_main >= 1 & df_work$LN_positive == 0) |
                       (df_work$SN_main == 0 & df_work$LN_positive >= 1), ]

cat("Patients with discordance between SN_main and LN_positive:\n")
if (nrow(discordant) > 0) {
  print(discordant[, c("Index", "SN_main", "LN_positive", "LN_total", "LNR")])
  cat("\n  Total discordant patients:", nrow(discordant), "\n")
  cat("  These patients are classified differently depending on the filter ",
      "strategy chosen in 00_config.R (NODE_POSITIVE_FILTER = '", 
      NODE_POSITIVE_FILTER, "').\n", sep = "")
} else {
  cat("  No discordances found. Both filter methods yield identical cohorts.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 9. SAVE OUTPUTS
# ═══════════════════════════════════════════════════════════════════════════════

section_header("SAVING OUTPUTS")

# --- Analytic dataset ---
save(df, df_all, df_raw, flow, ANALYTIC_TIER, 
     file = file.path(DIR_DATA, "analytic_dataset.RData"))
cat("Saved: analytic_dataset.RData\n")

write.csv(df, file.path(DIR_DATA, "analytic_dataset_nodepositive.csv"), 
          row.names = FALSE)
cat("Saved: analytic_dataset_nodepositive.csv\n")

# --- CONSORT flow summary ---
flow_df <- data.frame(
  Step = c("1. Total records in database",
           "2. Excluded: invalid LN (positive > total)",
           "3. Excluded: missing survival data",
           "4. Excluded: node-negative patients",
           "5. Excluded: missing/zero TLE",
           "Final analytic cohort"),
  N_remaining = c(flow$step1_total,
                  flow$step2_after_ln_valid,
                  flow$step3_after_surv,
                  flow$step4_node_positive,
                  flow$step5_final_analytic,
                  flow$step5_final_analytic),
  N_excluded = c(NA,
                 flow$excluded_invalid_ln,
                 flow$excluded_missing_surv,
                 flow$excluded_node_negative,
                 flow$excluded_tle_missing,
                 NA)
)
save_table(flow_df, "consort_flow")
print(flow_df)

# --- Missingness summary ---
save_table(missing_summary, "missingness_summary")

cat("\n", strrep("=", 78), "\n")
cat("  SCRIPT 01 COMPLETE\n")
cat("  Analytic cohort: N =", nrow(df), "| OS events:", n_os_events, 
    "| Tier:", ANALYTIC_TIER$label, "\n")
cat("  Filter method:", NODE_POSITIVE_FILTER, "\n")
cat(strrep("=", 78), "\n")
