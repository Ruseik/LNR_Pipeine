################################################################################
#  03_unadjusted_survival.R — Phase 2: Unadjusted Survival Analysis
#  Protocol Section 7.2 | SO1, SO2 foundational
#
#  Outputs: KM curves, log-rank tests, RMST, univariate Cox models,
#           LNR dose-response curve (martingale residual plot)
################################################################################

source(file.path("d:/Google_SSD_RAM/_OneDrive/Research",
                 "Project_ColoRectal/LNR_20260328/Analysis/00_config.R"))

load(file.path(DIR_DATA, "analytic_dataset.RData"))

section_header("SCRIPT 03: UNADJUSTED SURVIVAL ANALYSIS (Phase 2)")

cat("Analytic cohort: N =", nrow(df), "| OS events:", 
    sum(df$Outcome_OS == 1), "\n")
cat("Analytic tier:", ANALYTIC_TIER$label, "\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 1. KAPLAN-MEIER CURVES — OVERALL SURVIVAL BY N-STAGE
# ═══════════════════════════════════════════════════════════════════════════════

section_header("KM CURVES: OS BY N-STAGE")

# Fit KM
km_nstage <- survfit(Surv(Time_survival, Outcome_OS) ~ N_stage, data = df)
print(km_nstage)

# Log-rank test
lr_nstage <- survdiff(Surv(Time_survival, Outcome_OS) ~ N_stage, data = df)
cat("\nLog-rank test (N1 vs N2):\n")
print(lr_nstage)
lr_nstage_p <- 1 - pchisq(lr_nstage$chisq, df = length(lr_nstage$n) - 1)
cat("  p-value:", format_pval(lr_nstage_p), "\n")

# 1/3/5-year survival estimates
cat("\nSurvival estimates at landmark timepoints:\n")
surv_est_nstage <- summary(km_nstage, times = SURV_TIMEPOINTS_DAYS)
for (i in seq_along(SURV_TIMEPOINTS_YEARS)) {
  cat(sprintf("  %d-year survival:\n", SURV_TIMEPOINTS_YEARS[i]))
  idx <- which(surv_est_nstage$time == SURV_TIMEPOINTS_DAYS[i])
  if (length(idx) > 0) {
    for (j in idx) {
      cat(sprintf("    %s: %.1f%% (95%% CI: %.1f%% - %.1f%%)\n",
                  surv_est_nstage$strata[j],
                  surv_est_nstage$surv[j] * 100,
                  surv_est_nstage$lower[j] * 100,
                  surv_est_nstage$upper[j] * 100))
    }
  }
}

# KM plot
p_km_nstage <- ggsurvplot(
  km_nstage, data = df,
  palette = unname(COL_NSTAGE),
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.25,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = TRUE,
  conf.int.alpha = 0.15,
  xlab = "Time (days)",
  ylab = "Overall Survival Probability",
  title = "Overall Survival by N-Stage",
  subtitle = "Node-positive CRC | Log-rank test",
  legend.title = "N-Stage",
  legend.labs = c("N1", "N2"),
  ggtheme = theme_publication(),
  break.time.by = 365.25,
  surv.median.line = "hv",
  xlim = c(0, max(df$Time_survival, na.rm = TRUE))
)
pdf(file.path(DIR_FIGURES, "Fig_KM_OS_by_Nstage.pdf"), 
    width = FIG_WIDTH, height = FIG_HEIGHT + 1.5)
print(p_km_nstage)
dev.off()
png(file.path(DIR_FIGURES, "Fig_KM_OS_by_Nstage.png"),
    width = FIG_WIDTH, height = FIG_HEIGHT + 1.5, units = "in", res = FIG_DPI)
print(p_km_nstage)
dev.off()
cat("Saved: Fig_KM_OS_by_Nstage\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 2. KAPLAN-MEIER CURVES — OVERALL SURVIVAL BY LNR TIER
# ═══════════════════════════════════════════════════════════════════════════════

section_header("KM CURVES: OS BY LNR TIER")

km_lnr <- survfit(Surv(Time_survival, Outcome_OS) ~ LNR_tier, data = df)
print(km_lnr)

# Log-rank test
lr_lnr <- survdiff(Surv(Time_survival, Outcome_OS) ~ LNR_tier, data = df)
cat("\nLog-rank test (Low vs Intermediate vs High):\n")
print(lr_lnr)
lr_lnr_p <- 1 - pchisq(lr_lnr$chisq, df = length(lr_lnr$n) - 1)
cat("  p-value:", format_pval(lr_lnr_p), "\n")

# Pairwise log-rank with Bonferroni correction
cat("\nPairwise log-rank tests (Bonferroni corrected):\n")
if (nlevels(df$LNR_tier) > 1 && all(table(df$LNR_tier) > 0)) {
  pairwise_lr <- pairwise_survdiff(
    Surv(Time_survival, Outcome_OS) ~ LNR_tier,
    data = df, p.adjust.method = "bonferroni"
  )
  print(pairwise_lr)
}

# 1/3/5-year survival estimates
cat("\nSurvival estimates at landmark timepoints:\n")
surv_est_lnr <- summary(km_lnr, times = SURV_TIMEPOINTS_DAYS)
for (i in seq_along(SURV_TIMEPOINTS_YEARS)) {
  cat(sprintf("  %d-year survival:\n", SURV_TIMEPOINTS_YEARS[i]))
  idx <- which(surv_est_lnr$time == SURV_TIMEPOINTS_DAYS[i])
  if (length(idx) > 0) {
    for (j in idx) {
      cat(sprintf("    %s: %.1f%% (95%% CI: %.1f%% - %.1f%%)\n",
                  surv_est_lnr$strata[j],
                  surv_est_lnr$surv[j] * 100,
                  surv_est_lnr$lower[j] * 100,
                  surv_est_lnr$upper[j] * 100))
    }
  }
}

# KM plot
p_km_lnr <- ggsurvplot(
  km_lnr, data = df,
  palette = unname(COL_LNR_TIERS),
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.25,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = TRUE,
  conf.int.alpha = 0.15,
  xlab = "Time (days)",
  ylab = "Overall Survival Probability",
  title = "Overall Survival by LNR Tier (Rosenberg)",
  subtitle = "Node-positive CRC | Thresholds: 0.05 / 0.20",
  legend.title = "LNR Tier",
  legend.labs = LNR_TIER_LABELS,
  ggtheme = theme_publication(),
  break.time.by = 365.25,
  surv.median.line = "hv",
  xlim = c(0, max(df$Time_survival, na.rm = TRUE))
)
pdf(file.path(DIR_FIGURES, "Fig_KM_OS_by_LNR_tier.pdf"),
    width = FIG_WIDTH, height = FIG_HEIGHT + 1.5)
print(p_km_lnr)
dev.off()
png(file.path(DIR_FIGURES, "Fig_KM_OS_by_LNR_tier.png"),
    width = FIG_WIDTH, height = FIG_HEIGHT + 1.5, units = "in", res = FIG_DPI)
print(p_km_lnr)
dev.off()
cat("Saved: Fig_KM_OS_by_LNR_tier\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 3. RESTRICTED MEAN SURVIVAL TIME (RMST) AT 5 YEARS
# ═══════════════════════════════════════════════════════════════════════════════

section_header("RMST AT 5 YEARS")

# Determine truncation time — min of RMST_TAU_DAYS and max observed time
tau <- min(RMST_TAU_DAYS, max(df$Time_survival, na.rm = TRUE))
cat("RMST truncation time (tau):", round(tau, 0), "days (", 
    round(tau / 365.25, 2), "years)\n\n")

# RMST by N-stage
cat("RMST by N-Stage:\n")
nstage_levels <- levels(df$N_stage)
for (lev in nstage_levels) {
  sub <- df[df$N_stage == lev, ]
  if (nrow(sub) > 0 && sum(sub$Outcome_OS == 1) > 0) {
    km_sub <- survfit(Surv(Time_survival, Outcome_OS) ~ 1, data = sub)
    # Compute RMST as area under KM curve up to tau
    rmst_val <- summary(km_sub, rmean = tau)$table["*rmean"]
    rmst_se <- summary(km_sub, rmean = tau)$table["*se(rmean)"]
    cat(sprintf("  %s: RMST = %.0f days (%.1f years) [SE: %.0f]\n",
                lev, rmst_val, rmst_val / 365.25, rmst_se))
  }
}

# RMST by LNR tier
cat("\nRMST by LNR Tier:\n")
lnr_levels <- levels(df$LNR_tier)
for (lev in lnr_levels) {
  sub <- df[df$LNR_tier == lev, ]
  if (nrow(sub) > 0 && sum(sub$Outcome_OS == 1) > 0) {
    km_sub <- survfit(Surv(Time_survival, Outcome_OS) ~ 1, data = sub)
    rmst_val <- summary(km_sub, rmean = tau)$table["*rmean"]
    rmst_se <- summary(km_sub, rmean = tau)$table["*se(rmean)"]
    cat(sprintf("  %s: RMST = %.0f days (%.1f years) [SE: %.0f]\n",
                lev, rmst_val, rmst_val / 365.25, rmst_se))
  } else {
    cat(sprintf("  %s: insufficient data (n = %d, events = %d)\n",
                lev, nrow(sub), sum(sub$Outcome_OS == 1)))
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# 4. CONTINUOUS LNR DOSE-RESPONSE (Martingale Residual Plot)
# ═══════════════════════════════════════════════════════════════════════════════

section_header("LNR DOSE-RESPONSE ASSESSMENT")

# Martingale residuals from null Cox model
null_cox <- coxph(Surv(Time_survival, Outcome_OS) ~ 1, data = df)
df$mart_resid <- residuals(null_cox, type = "martingale")

p_dose_response <- ggplot(df, aes(x = LNR_prop, y = mart_resid)) +
  geom_point(alpha = 0.5, colour = "grey40", size = 2) +
  geom_smooth(method = "loess", se = TRUE, colour = "#B2182B", 
              fill = "#B2182B", alpha = 0.15, linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = c(LNR_CUT_LOW_INTER, LNR_CUT_INTER_HIGH),
             linetype = "dotted", colour = c("#4DAF4A", "#FF7F00"),
             linewidth = 0.6) +
  labs(title = "LNR Dose-Response: Martingale Residuals from Null Cox Model",
       subtitle = "LOESS smoothed | Rosenberg thresholds indicated",
       x = "Lymph Node Ratio (proportion)",
       y = "Martingale Residual") +
  theme_publication()
save_figure(p_dose_response, "Fig_LNR_dose_response_martingale")

cat("Visual assessment: Check whether the LOESS curve shows a\n")
cat("monotonically increasing relationship (linear dose-response)\n")
cat("or a non-linear (threshold/plateau) pattern.\n")
cat("This determines whether RCS parametrisation is justified.\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 5. UNIVARIATE COX MODELS
# ═══════════════════════════════════════════════════════════════════════════════

section_header("UNIVARIATE COX MODELS")

# --- 5.1 N-Stage (layer 1, Model A) ---
cat("--- Model A: N-Stage (univariate) ---\n")
cox_nstage_uni <- coxph(Surv(Time_survival, Outcome_OS) ~ N_stage, data = df)
print(summary(cox_nstage_uni))

# --- 5.2 LNR continuous (linear) ---
cat("\n--- LNR continuous (linear, univariate) ---\n")
cox_lnr_linear <- coxph(Surv(Time_survival, Outcome_OS) ~ LNR_prop, data = df)
print(summary(cox_lnr_linear))

# --- 5.3 LNR continuous with RCS (Layer 1, Model B) ---
# Only if tier permits RCS
if (ANALYTIC_TIER$permit_rcs || ANALYTIC_TIER$tier_code >= 2) {
  cat("\n--- Model B: LNR continuous (RCS, 3 knots, univariate) ---\n")
  # Set data distribution for rms
  dd <- datadist(df)
  options(datadist = "dd")
  
  cox_lnr_rcs <- cph(Surv(Time_survival, Outcome_OS) ~ rcs(LNR_prop, 3),
                      data = df, x = TRUE, y = TRUE, surv = TRUE)
  print(cox_lnr_rcs)
  cat("\nAnova (non-linearity test):\n")
  print(anova(cox_lnr_rcs))
  
  # RCS association plot (log-HR vs LNR)
  p_rcs <- ggplot(Predict(cox_lnr_rcs, LNR_prop, fun = exp),
                  aes(x = LNR_prop, y = yhat)) +
    geom_line(colour = "#B2182B", linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, fill = "#B2182B") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    labs(title = "LNR–Mortality Association (RCS, 3 knots)",
         subtitle = "Hazard ratio relative to reference LNR",
         x = "Lymph Node Ratio (proportion)",
         y = "Hazard Ratio") +
    theme_publication()
  save_figure(p_rcs, "Fig_LNR_RCS_association_plot")
} else {
  cat("\n  RCS not permitted at this analytic tier. Using linear LNR.\n")
  # For consistency, define Model B as linear LNR
  dd <- datadist(df)
  options(datadist = "dd")
  cox_lnr_rcs <- cph(Surv(Time_survival, Outcome_OS) ~ LNR_prop,
                      data = df, x = TRUE, y = TRUE, surv = TRUE)
  print(cox_lnr_rcs)
}

# --- 5.4 LNR categorical tiers ---
cat("\n--- LNR categorical tiers (Rosenberg, univariate) ---\n")
cox_lnr_cat <- coxph(Surv(Time_survival, Outcome_OS) ~ LNR_tier, data = df)
print(summary(cox_lnr_cat))

# --- 5.5 All other covariates (univariate) ---
cat("\n--- Other covariates (univariate) ---\n")
uni_vars <- c("AG", "Sex_f", "T_stage_group", "Grade_f", "LVI_f", 
              "LN_total", "TLE_adequate")

uni_results <- data.frame(
  Variable = character(),
  HR = numeric(),
  HR_lower = numeric(),
  HR_upper = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (v in uni_vars) {
  formula_v <- as.formula(paste("Surv(Time_survival, Outcome_OS) ~", v))
  fit_v <- tryCatch(
    coxph(formula_v, data = df),
    error = function(e) NULL
  )
  if (!is.null(fit_v)) {
    s <- summary(fit_v)
    for (r in seq_len(nrow(s$coefficients))) {
      uni_results <- rbind(uni_results, data.frame(
        Variable = rownames(s$coefficients)[r],
        HR = round(s$conf.int[r, "exp(coef)"], 3),
        HR_lower = round(s$conf.int[r, "lower .95"], 3),
        HR_upper = round(s$conf.int[r, "upper .95"], 3),
        p_value = round(s$coefficients[r, "Pr(>|z|)"], 4),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Add LNR and N-stage
s_ns <- summary(cox_nstage_uni)
for (r in seq_len(nrow(s_ns$coefficients))) {
  uni_results <- rbind(uni_results, data.frame(
    Variable = rownames(s_ns$coefficients)[r],
    HR = round(s_ns$conf.int[r, "exp(coef)"], 3),
    HR_lower = round(s_ns$conf.int[r, "lower .95"], 3),
    HR_upper = round(s_ns$conf.int[r, "upper .95"], 3),
    p_value = round(s_ns$coefficients[r, "Pr(>|z|)"], 4),
    stringsAsFactors = FALSE
  ))
}

s_lnr <- summary(cox_lnr_linear)
uni_results <- rbind(uni_results, data.frame(
  Variable = "LNR_prop (continuous)",
  HR = round(s_lnr$conf.int[1, "exp(coef)"], 3),
  HR_lower = round(s_lnr$conf.int[1, "lower .95"], 3),
  HR_upper = round(s_lnr$conf.int[1, "upper .95"], 3),
  p_value = round(s_lnr$coefficients[1, "Pr(>|z|)"], 4),
  stringsAsFactors = FALSE
))

cat("\nUnivariate Cox regression — all covariates:\n")
print(uni_results)
save_table(uni_results, "univariate_cox_results")

# ═══════════════════════════════════════════════════════════════════════════════
# 6. DISEASE-SPECIFIC SURVIVAL — KM (if applicable)
# ═══════════════════════════════════════════════════════════════════════════════

section_header("KM CURVES: DSS (DISEASE-SPECIFIC SURVIVAL)")

n_dss_events <- sum(df$Outcome_DSS == 1, na.rm = TRUE)
cat("DSS events:", n_dss_events, "\n")

if (n_dss_events >= 5) {
  # KM by N-stage
  km_dss_nstage <- survfit(Surv(Time_survival, Outcome_DSS) ~ N_stage, data = df)
  p_km_dss_ns <- ggsurvplot(
    km_dss_nstage, data = df, palette = unname(COL_NSTAGE),
    risk.table = TRUE, risk.table.height = 0.25,
    pval = TRUE, pval.method = TRUE, conf.int = TRUE, conf.int.alpha = 0.15,
    xlab = "Time (days)", ylab = "Disease-Specific Survival Probability",
    title = "Disease-Specific Survival by N-Stage",
    legend.title = "N-Stage", legend.labs = c("N1", "N2"),
    ggtheme = theme_publication(), break.time.by = 365.25
  )
  pdf(file.path(DIR_FIGURES, "Fig_KM_DSS_by_Nstage.pdf"),
      width = FIG_WIDTH, height = FIG_HEIGHT + 1.5)
  print(p_km_dss_ns)
  dev.off()
  png(file.path(DIR_FIGURES, "Fig_KM_DSS_by_Nstage.png"),
      width = FIG_WIDTH, height = FIG_HEIGHT + 1.5, units = "in", res = FIG_DPI)
  print(p_km_dss_ns)
  dev.off()
  
  # KM by LNR tier
  km_dss_lnr <- survfit(Surv(Time_survival, Outcome_DSS) ~ LNR_tier, data = df)
  p_km_dss_lnr <- ggsurvplot(
    km_dss_lnr, data = df, palette = unname(COL_LNR_TIERS),
    risk.table = TRUE, risk.table.height = 0.25,
    pval = TRUE, pval.method = TRUE, conf.int = TRUE, conf.int.alpha = 0.15,
    xlab = "Time (days)", ylab = "Disease-Specific Survival Probability",
    title = "Disease-Specific Survival by LNR Tier",
    legend.title = "LNR Tier", legend.labs = LNR_TIER_LABELS,
    ggtheme = theme_publication(), break.time.by = 365.25
  )
  pdf(file.path(DIR_FIGURES, "Fig_KM_DSS_by_LNR_tier.pdf"),
      width = FIG_WIDTH, height = FIG_HEIGHT + 1.5)
  print(p_km_dss_lnr)
  dev.off()
  png(file.path(DIR_FIGURES, "Fig_KM_DSS_by_LNR_tier.png"),
      width = FIG_WIDTH, height = FIG_HEIGHT + 1.5, units = "in", res = FIG_DPI)
  print(p_km_dss_lnr)
  dev.off()
  
  cat("Saved: Fig_KM_DSS_by_Nstage, Fig_KM_DSS_by_LNR_tier\n")
} else {
  cat("  Insufficient DSS events for KM analysis.\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 7. SAVE SESSION OBJECTS
# ═══════════════════════════════════════════════════════════════════════════════

save(km_nstage, km_lnr, cox_nstage_uni, cox_lnr_linear, cox_lnr_rcs,
     cox_lnr_cat, uni_results, lr_nstage_p, lr_lnr_p,
     file = file.path(DIR_MODELS, "phase2_unadjusted_results.RData"))

cat("\n", strrep("=", 78), "\n")
cat("  SCRIPT 03 COMPLETE\n")
cat(strrep("=", 78), "\n")
