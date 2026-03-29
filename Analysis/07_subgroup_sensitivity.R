################################################################################
#  07_subgroup_sensitivity.R — Phase 5: Subgroup & Sensitivity Analyses
#  Protocol Section 7.6 | SO5
#
#  Pre-specified subgroups: TLE <12 vs >=12, T3/T4 vs T1/T2, N1 vs N2
#  Pre-specified sensitivities: landmark, binary TLE, alt LNR cutpoints
################################################################################

source(file.path("d:/Google_SSD_RAM/_OneDrive/Research",
                 "Project_ColoRectal/LNR_20260328/Analysis/00_config.R"))

load(file.path(DIR_DATA, "analytic_dataset.RData"))
load(file.path(DIR_MODELS, "phase3_cox_models.RData"))

section_header("SCRIPT 07: SUBGROUP & SENSITIVITY ANALYSES (Phase 5)")

cat("Analytic cohort: N =", nrow(df), "| OS events:", 
    sum(df$Outcome_OS == 1), "\n")
cat("Analytic tier:", ANALYTIC_TIER$label, "\n\n")

# Re-establish rms datadist
dd <- datadist(df)
options(datadist = "dd")

# Helper to compute C-stat for a model within a subgroup
compute_c_subgroup <- function(formula, data, label) {
  tryCatch({
    fit <- coxph(formula, data = data)
    lp <- predict(fit, type = "lp")
    conc <- concordance(Surv(Time_survival, Outcome_OS) ~ lp, 
                        data = data, timewt = "n/G2")
    data.frame(
      Subgroup = label,
      N = nrow(data),
      Events = sum(data$Outcome_OS == 1),
      C_statistic = round(conc$concordance, 4),
      SE = round(sqrt(conc$var), 4),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      Subgroup = label, N = nrow(data),
      Events = sum(data$Outcome_OS == 1),
      C_statistic = NA, SE = NA, stringsAsFactors = FALSE
    )
  })
}

# ═══════════════════════════════════════════════════════════════════════════════
# 1. SUBGROUP ANALYSIS: TLE < 12 vs >= 12 (H4 — Primary subgroup)
# ═══════════════════════════════════════════════════════════════════════════════

section_header("SUBGROUP: TLE ADEQUACY (< 12 vs >= 12)")

cat("This is the primary subgroup of mechanistic interest (H4).\n")
cat("Tests whether LNR advantage is amplified in sub-optimal harvest.\n\n")

# Define formulas based on tier
if (ANALYTIC_TIER$permit_rcs) {
  formula_nstage <- Surv(Time_survival, Outcome_OS) ~ N_stage
  formula_lnr <- Surv(Time_survival, Outcome_OS) ~ LNR_prop
} else {
  formula_nstage <- Surv(Time_survival, Outcome_OS) ~ N_stage
  formula_lnr <- Surv(Time_survival, Outcome_OS) ~ LNR_prop
}

subgroup_results_tle <- data.frame()

for (tle_level in c("Inadequate (<12)", "Adequate (>=12)")) {
  sub_data <- df[df$TLE_adequate == tle_level, ]
  cat(sprintf("\n━━━ TLE %s (N = %d, events = %d) ━━━\n",
              tle_level, nrow(sub_data), sum(sub_data$Outcome_OS == 1)))
  
  if (nrow(sub_data) >= 10 && sum(sub_data$Outcome_OS == 1) >= 3) {
    # N-stage model C-stat
    c_ns <- compute_c_subgroup(formula_nstage, sub_data, 
                                paste("N-stage |", tle_level))
    # LNR model C-stat
    c_lnr <- compute_c_subgroup(formula_lnr, sub_data, 
                                 paste("LNR |", tle_level))
    
    subgroup_results_tle <- rbind(subgroup_results_tle, c_ns, c_lnr)
    
    delta_c <- c_lnr$C_statistic - c_ns$C_statistic
    cat(sprintf("  N-stage C = %.4f | LNR C = %.4f | ΔC = %.4f\n",
                c_ns$C_statistic, c_lnr$C_statistic, delta_c))
    
    # KM curves within subgroup
    if (sum(sub_data$Outcome_OS == 1) >= 3 && 
        nlevels(droplevels(sub_data$LNR_tier)) >= 2) {
      km_sub <- survfit(Surv(Time_survival, Outcome_OS) ~ LNR_tier, data = sub_data)
      p_km_sub <- ggsurvplot(
        km_sub, data = sub_data,
        palette = unname(COL_LNR_TIERS),
        pval = TRUE, conf.int = FALSE,
        xlab = "Time (days)", ylab = "Overall Survival",
        title = paste("OS by LNR Tier — TLE", tle_level),
        ggtheme = theme_publication(), break.time.by = 365.25,
        risk.table = TRUE, risk.table.height = 0.25
      )
      pdf(file.path(DIR_FIGURES, 
                     paste0("Fig_KM_LNR_TLE_", gsub("[^a-zA-Z0-9]", "", tle_level), ".pdf")),
          width = FIG_WIDTH, height = FIG_HEIGHT + 1.5)
      print(p_km_sub)
      dev.off()
    }
  } else {
    cat("  Insufficient data for subgroup analysis.\n")
  }
}

cat("\nSubgroup comparison summary (TLE adequacy):\n")
print(subgroup_results_tle)
save_table(subgroup_results_tle, "subgroup_TLE_adequacy")

# Formal interaction test (only at tier 4)
if (ANALYTIC_TIER$permit_interaction) {
  cat("\n━━━ Formal Interaction Test: LNR × TLE adequacy ━━━\n")
  int_model <- coxph(Surv(Time_survival, Outcome_OS) ~ 
                       LNR_prop * TLE_adequate_binary, data = df)
  print(summary(int_model))
  int_anova <- anova(int_model)
  print(int_anova)
} else {
  cat("\nFormal interaction test NOT permitted (requires >= 200 events).\n")
  cat("Subgroup results above are DESCRIPTIVE and EXPLORATORY only.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 2. SUBGROUP: T-STAGE (T3/T4 vs T1/T2)
# ═══════════════════════════════════════════════════════════════════════════════

section_header("SUBGROUP: T-STAGE (T3/T4 vs T1/T2)")

df$T_high <- ifelse(df$ST_main >= 3, "T3/T4", "T1/T2")

subgroup_results_tstage <- data.frame()
for (t_level in c("T1/T2", "T3/T4")) {
  sub_data <- df[df$T_high == t_level, ]
  cat(sprintf("\n━━━ %s (N = %d, events = %d) ━━━\n",
              t_level, nrow(sub_data), sum(sub_data$Outcome_OS == 1)))
  
  if (nrow(sub_data) >= 10 && sum(sub_data$Outcome_OS == 1) >= 3) {
    c_ns <- compute_c_subgroup(formula_nstage, sub_data, paste("N-stage |", t_level))
    c_lnr <- compute_c_subgroup(formula_lnr, sub_data, paste("LNR |", t_level))
    subgroup_results_tstage <- rbind(subgroup_results_tstage, c_ns, c_lnr)
    delta_c <- c_lnr$C_statistic - c_ns$C_statistic
    cat(sprintf("  N-stage C = %.4f | LNR C = %.4f | ΔC = %.4f\n",
                c_ns$C_statistic, c_lnr$C_statistic, delta_c))
  } else {
    cat("  Insufficient data.\n")
  }
}

print(subgroup_results_tstage)
save_table(subgroup_results_tstage, "subgroup_T_stage")

# ═══════════════════════════════════════════════════════════════════════════════
# 3. SUBGROUP: N1 vs N2
# ═══════════════════════════════════════════════════════════════════════════════

section_header("SUBGROUP: N1 vs N2")

subgroup_results_nstage <- data.frame()
for (n_level in c("N1", "N2")) {
  sub_data <- df[df$N_stage == n_level, ]
  cat(sprintf("\n━━━ %s (N = %d, events = %d) ━━━\n",
              n_level, nrow(sub_data), sum(sub_data$Outcome_OS == 1)))
  
  if (nrow(sub_data) >= 10 && sum(sub_data$Outcome_OS == 1) >= 3) {
    c_lnr <- compute_c_subgroup(formula_lnr, sub_data, paste("LNR |", n_level))
    subgroup_results_nstage <- rbind(subgroup_results_nstage, c_lnr)
    cat(sprintf("  LNR C-statistic within %s stratum: %.4f\n",
                n_level, c_lnr$C_statistic))
  } else {
    cat("  Insufficient data.\n")
  }
}

print(subgroup_results_nstage)
save_table(subgroup_results_nstage, "subgroup_N_stage")

# ═══════════════════════════════════════════════════════════════════════════════
# 4. SENSITIVITY ANALYSIS: LANDMARK AT 6 MONTHS
# ═══════════════════════════════════════════════════════════════════════════════

section_header("SENSITIVITY: 6-MONTH LANDMARK ANALYSIS")

cat("Excluding events and observations within first", LANDMARK_DAYS, "days.\n")
cat("Purpose: Remove early surgical mortality from prognostic signal.\n\n")

df_landmark <- df[df$Time_survival > LANDMARK_DAYS, ]
df_landmark$Time_landmark <- df_landmark$Time_survival - LANDMARK_DAYS

cat("After landmark exclusion: N =", nrow(df_landmark), 
    "| Events:", sum(df_landmark$Outcome_OS == 1), "\n")
cat("Excluded:", nrow(df) - nrow(df_landmark), "patients with", 
    sum(df$Outcome_OS == 1 & df$Time_survival <= LANDMARK_DAYS), 
    "early deaths.\n\n")

if (nrow(df_landmark) >= 20 && sum(df_landmark$Outcome_OS == 1) >= 5) {
  # Refit univariate models on landmark cohort
  cox_ns_lm <- coxph(Surv(Time_landmark, Outcome_OS) ~ N_stage, data = df_landmark)
  cox_lnr_lm <- coxph(Surv(Time_landmark, Outcome_OS) ~ LNR_prop, data = df_landmark)
  
  cat("N-stage (landmark):\n")
  print(summary(cox_ns_lm))
  cat("\nLNR continuous (landmark):\n")
  print(summary(cox_lnr_lm))
  
  # C-statistics
  c_ns_lm <- compute_c_subgroup(
    Surv(Time_landmark, Outcome_OS) ~ N_stage, df_landmark, "N-stage (landmark)")
  c_lnr_lm <- compute_c_subgroup(
    Surv(Time_landmark, Outcome_OS) ~ LNR_prop, df_landmark, "LNR (landmark)")
  
  landmark_results <- rbind(c_ns_lm, c_lnr_lm)
  print(landmark_results)
  save_table(landmark_results, "sensitivity_landmark_6mo")
  
  # KM by LNR tier (landmark)
  km_lnr_lm <- survfit(Surv(Time_landmark, Outcome_OS) ~ LNR_tier, 
                         data = df_landmark)
  p_km_lm <- ggsurvplot(
    km_lnr_lm, data = df_landmark, palette = unname(COL_LNR_TIERS),
    pval = TRUE, conf.int = FALSE, risk.table = TRUE, risk.table.height = 0.25,
    xlab = "Time from Landmark (days)", ylab = "Overall Survival",
    title = "OS by LNR Tier — 6-Month Landmark Analysis",
    ggtheme = theme_publication(), break.time.by = 365.25
  )
  pdf(file.path(DIR_FIGURES, "Fig_KM_landmark_6mo.pdf"),
      width = FIG_WIDTH, height = FIG_HEIGHT + 1.5)
  print(p_km_lm)
  dev.off()
} else {
  cat("  Insufficient data after landmark exclusion.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 5. SENSITIVITY: BINARY TLE IN MVA
# ═══════════════════════════════════════════════════════════════════════════════

if (ANALYTIC_TIER$permit_mva) {
  section_header("SENSITIVITY: BINARY TLE (< 12 vs >= 12) IN MVA")
  
  cat("Replacing continuous/RCS TLE with binary (< 12 vs >= 12).\n\n")
  
  # Model C with binary TLE
  cox_C_binTLE <- coxph(Surv(Time_survival, Outcome_OS) ~ 
                          N_stage + TLE_adequate_binary + AG + Sex_f + 
                          T_stage_group + Grade_f + LVI_f,
                        data = df)
  cat("Model C (binary TLE):\n")
  print(summary(cox_C_binTLE))
  
  # Model D with binary TLE
  cox_D_binTLE <- coxph(Surv(Time_survival, Outcome_OS) ~ 
                          LNR_prop + TLE_adequate_binary + AG + Sex_f + 
                          T_stage_group + Grade_f + LVI_f,
                        data = df)
  cat("\nModel D (binary TLE):\n")
  print(summary(cox_D_binTLE))
  
  # C-statistics
  c_C_bin <- compute_c_subgroup(
    Surv(Time_survival, Outcome_OS) ~ N_stage + TLE_adequate_binary + AG + Sex_f + 
      T_stage_group + Grade_f + LVI_f,
    df, "Model C (binary TLE)")
  c_D_bin <- compute_c_subgroup(
    Surv(Time_survival, Outcome_OS) ~ LNR_prop + TLE_adequate_binary + AG + Sex_f + 
      T_stage_group + Grade_f + LVI_f,
    df, "Model D (binary TLE)")
  
  binTLE_results <- rbind(c_C_bin, c_D_bin)
  binTLE_results$Delta_C <- c(NA, c_D_bin$C_statistic - c_C_bin$C_statistic)
  print(binTLE_results)
  save_table(binTLE_results, "sensitivity_binary_TLE")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 6. SENSITIVITY: ALTERNATIVE LNR CUTPOINTS
# ═══════════════════════════════════════════════════════════════════════════════

section_header("SENSITIVITY: ALTERNATIVE LNR CUTPOINTS")

# --- 6.1 Cutpoints 0.05/0.40 ---
cat("━━━ LNR tiers with 0.05/0.40 cutpoints ━━━\n")
km_alt1 <- survfit(Surv(Time_survival, Outcome_OS) ~ LNR_tier_alt1, data = df)
cat("Distribution:\n")
print(table(df$LNR_tier_alt1))
lr_alt1 <- survdiff(Surv(Time_survival, Outcome_OS) ~ LNR_tier_alt1, data = df)
cat("Log-rank p:", format_pval(1 - pchisq(lr_alt1$chisq, 
                                           df = length(lr_alt1$n) - 1)), "\n\n")

# --- 6.2 Binary cutpoint at 0.20 ---
cat("━━━ LNR binary with 0.20 cutpoint ━━━\n")
km_alt2 <- survfit(Surv(Time_survival, Outcome_OS) ~ LNR_tier_binary, data = df)
cat("Distribution:\n")
print(table(df$LNR_tier_binary))
lr_alt2 <- survdiff(Surv(Time_survival, Outcome_OS) ~ LNR_tier_binary, data = df)
cat("Log-rank p:", format_pval(1 - pchisq(lr_alt2$chisq, 
                                           df = length(lr_alt2$n) - 1)), "\n")

# KM plot for binary threshold
p_km_binary <- ggsurvplot(
  km_alt2, data = df,
  palette = c("#2166AC", "#B2182B"),
  risk.table = TRUE, risk.table.height = 0.25,
  pval = TRUE, pval.method = TRUE,
  conf.int = TRUE, conf.int.alpha = 0.15,
  xlab = "Time (days)", ylab = "Overall Survival Probability",
  title = "OS by LNR Binary (0.20 Threshold) — Sensitivity",
  legend.title = "LNR Group",
  ggtheme = theme_publication(), break.time.by = 365.25
)
pdf(file.path(DIR_APPENDIX, "Fig_KM_LNR_binary_020.pdf"),
    width = FIG_WIDTH, height = FIG_HEIGHT + 1.5)
print(p_km_binary)
dev.off()

# --- 6.3 Compare C-statistics across LNR categorisations ---
cat("\n━━━ C-statistic comparison across LNR categorisations ━━━\n")
lnr_cuts_compare <- data.frame()
for (tier_var in c("LNR_tier", "LNR_tier_alt1", "LNR_tier_binary")) {
  fml <- as.formula(paste("Surv(Time_survival, Outcome_OS) ~", tier_var))
  c_res <- compute_c_subgroup(fml, df, tier_var)
  lnr_cuts_compare <- rbind(lnr_cuts_compare, c_res)
}
# Add continuous
c_cont <- compute_c_subgroup(
  Surv(Time_survival, Outcome_OS) ~ LNR_prop, df, "LNR_prop (continuous)")
lnr_cuts_compare <- rbind(lnr_cuts_compare, c_cont)

print(lnr_cuts_compare)
save_table(lnr_cuts_compare, "sensitivity_LNR_cutpoints_comparison")

# ═══════════════════════════════════════════════════════════════════════════════
# 7. COMPILE ALL SUBGROUP/SENSITIVITY RESULTS
# ═══════════════════════════════════════════════════════════════════════════════

section_header("SUMMARY TABLE: ALL SUBGROUP & SENSITIVITY ANALYSES")

all_subgroup <- rbind(
  if (nrow(subgroup_results_tle) > 0) 
    cbind(Analysis = "TLE adequacy subgroup", subgroup_results_tle) else NULL,
  if (nrow(subgroup_results_tstage) > 0) 
    cbind(Analysis = "T-stage subgroup", subgroup_results_tstage) else NULL,
  if (nrow(subgroup_results_nstage) > 0) 
    cbind(Analysis = "N-stage subgroup", subgroup_results_nstage) else NULL,
  if (exists("landmark_results")) 
    cbind(Analysis = "6-month landmark", landmark_results) else NULL
)

if (!is.null(all_subgroup) && nrow(all_subgroup) > 0) {
  cat("Complete subgroup & sensitivity results:\n")
  print(all_subgroup)
  save_table(all_subgroup, "all_subgroup_sensitivity_results")
}

save(subgroup_results_tle, subgroup_results_tstage, subgroup_results_nstage,
     lnr_cuts_compare,
     file = file.path(DIR_MODELS, "phase5_subgroup_sensitivity.RData"))

cat("\n", strrep("=", 78), "\n")
cat("  SCRIPT 07 COMPLETE\n")
cat(strrep("=", 78), "\n")
