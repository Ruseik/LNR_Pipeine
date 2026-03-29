################################################################################
#  04_cox_models.R — Phase 3: Three-Layer Cox Model Comparison
#  Protocol Section 7.3 | SO2, SO3
#
#  Implements the six-model symmetric comparison framework:
#    Layer 1 (Unadjusted): Model A (N-stage) vs Model B (LNR RCS)
#    Layer 2 (TLE-adjusted): Model C2 vs Model D2
#    Layer 3 (Full MVA):     Model C vs Model D  (PRIMARY)
#
#  Model complexity adapts automatically to the event-count gating tier.
################################################################################

source(file.path("d:/Google_SSD_RAM/_OneDrive/Research",
                 "Project_ColoRectal/LNR_20260328/Analysis/00_config.R"))

load(file.path(DIR_DATA, "analytic_dataset.RData"))

section_header("SCRIPT 04: THREE-LAYER COX MODEL COMPARISON (Phase 3)")

cat("Analytic cohort: N =", nrow(df), "| OS events:", 
    sum(df$Outcome_OS == 1), "\n")
cat("Analytic tier:", ANALYTIC_TIER$label, "\n")
cat("Permit MVA:", ANALYTIC_TIER$permit_mva, 
    "| Permit RCS:", ANALYTIC_TIER$permit_rcs, "\n\n")

# Set up rms data distribution
dd <- datadist(df)
options(datadist = "dd")

# Store all models in a named list
models <- list()
model_summary <- data.frame(
  Model = character(),
  Layer = character(),
  Description = character(),
  df_used = integer(),
  AIC = numeric(),
  BIC = numeric(),
  LogLik = numeric(),
  stringsAsFactors = FALSE
)

# ═══════════════════════════════════════════════════════════════════════════════
# 1. LAYER 1 — UNADJUSTED MODELS
# ═══════════════════════════════════════════════════════════════════════════════

section_header("LAYER 1: UNADJUSTED MODELS")

# --- Model A: N-stage only ---
cat("━━━ Model A: N-stage only ━━━\n")
models$A <- cph(Surv(Time_survival, Outcome_OS) ~ N_stage,
                data = df, x = TRUE, y = TRUE, surv = TRUE)
print(models$A)
cat("\n")

model_summary <- rbind(model_summary, data.frame(
  Model = "A", Layer = "1 (Unadjusted)",
  Description = "N-stage only",
  df_used = models$A$stats["d.f."],
  AIC = AIC(models$A), BIC = BIC(models$A),
  LogLik = models$A$loglik[2]
))

# --- Model B: LNR (continuous) ---
cat("━━━ Model B: LNR continuous ━━━\n")
if (ANALYTIC_TIER$permit_rcs) {
  models$B <- cph(Surv(Time_survival, Outcome_OS) ~ rcs(LNR_prop, 3),
                  data = df, x = TRUE, y = TRUE, surv = TRUE)
  cat("  (RCS with 3 knots)\n")
} else {
  models$B <- cph(Surv(Time_survival, Outcome_OS) ~ LNR_prop,
                  data = df, x = TRUE, y = TRUE, surv = TRUE)
  cat("  (Linear — RCS not permitted at this tier)\n")
}
print(models$B)
cat("\nAnova:\n")
print(anova(models$B))

model_summary <- rbind(model_summary, data.frame(
  Model = "B", Layer = "1 (Unadjusted)",
  Description = ifelse(ANALYTIC_TIER$permit_rcs, "LNR RCS(3)", "LNR linear"),
  df_used = models$B$stats["d.f."],
  AIC = AIC(models$B), BIC = BIC(models$B),
  LogLik = models$B$loglik[2]
))

# ═══════════════════════════════════════════════════════════════════════════════
# 2. LAYER 2 — TLE-ADJUSTED MODELS (if MVA permitted)
# ═══════════════════════════════════════════════════════════════════════════════

if (ANALYTIC_TIER$permit_mva) {
  section_header("LAYER 2: TLE-ADJUSTED MODELS")
  
  # --- Model C2: N-stage + TLE ---
  cat("━━━ Model C2: N-stage + TLE ━━━\n")
  if (ANALYTIC_TIER$permit_rcs) {
    models$C2 <- cph(Surv(Time_survival, Outcome_OS) ~ N_stage + rcs(LN_total, 3),
                     data = df, x = TRUE, y = TRUE, surv = TRUE)
    cat("  (TLE modelled with RCS, 3 knots)\n")
  } else {
    models$C2 <- cph(Surv(Time_survival, Outcome_OS) ~ N_stage + LN_total,
                     data = df, x = TRUE, y = TRUE, surv = TRUE)
    cat("  (TLE modelled linearly)\n")
  }
  print(models$C2)
  cat("\nAnova:\n")
  print(anova(models$C2))
  
  model_summary <- rbind(model_summary, data.frame(
    Model = "C2", Layer = "2 (TLE-adjusted)",
    Description = ifelse(ANALYTIC_TIER$permit_rcs, 
                         "N-stage + TLE RCS(3)", "N-stage + TLE linear"),
    df_used = models$C2$stats["d.f."],
    AIC = AIC(models$C2), BIC = BIC(models$C2),
    LogLik = models$C2$loglik[2]
  ))
  
  # --- Model D2: LNR + TLE ---
  cat("\n━━━ Model D2: LNR + TLE ━━━\n")
  if (ANALYTIC_TIER$permit_rcs) {
    models$D2 <- cph(Surv(Time_survival, Outcome_OS) ~ 
                       rcs(LNR_prop, 3) + rcs(LN_total, 3),
                     data = df, x = TRUE, y = TRUE, surv = TRUE)
    cat("  (Both LNR and TLE modelled with RCS, 3 knots each)\n")
  } else {
    models$D2 <- cph(Surv(Time_survival, Outcome_OS) ~ LNR_prop + LN_total,
                     data = df, x = TRUE, y = TRUE, surv = TRUE)
    cat("  (Both LNR and TLE modelled linearly)\n")
  }
  print(models$D2)
  cat("\nAnova:\n")
  print(anova(models$D2))
  
  model_summary <- rbind(model_summary, data.frame(
    Model = "D2", Layer = "2 (TLE-adjusted)",
    Description = ifelse(ANALYTIC_TIER$permit_rcs, 
                         "LNR RCS(3) + TLE RCS(3)", "LNR linear + TLE linear"),
    df_used = models$D2$stats["d.f."],
    AIC = AIC(models$D2), BIC = BIC(models$D2),
    LogLik = models$D2$loglik[2]
  ))
  
} else {
  cat("\n  Layer 2 SKIPPED — MVA not permitted at Tier", 
      ANALYTIC_TIER$tier_code, "\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 3. LAYER 3 — FULL MVA MODELS (PRIMARY COMPARISON)
# ═══════════════════════════════════════════════════════════════════════════════

if (ANALYTIC_TIER$permit_mva) {
  section_header("LAYER 3: FULL MVA MODELS (PRIMARY)")
  
  # Define covariate block — forced entry, no stepwise selection (per protocol)
  # Available covariates: AG, Sex_f, T_stage_group, Grade_f, LVI_f, LN_total
  # Missing from protocol: PNI, resection margin, tumour site, calendar period
  
  # --- Model C: N-stage + TLE + all covariates (Primary N-stage model) ---
  cat("━━━ Model C: N-stage Full MVA (PRIMARY) ━━━\n")
  cat("  Covariates: N-stage + TLE + Age + Sex + T-stage + Grade + LVI\n")
  cat("  (PNI, margin, site, calendar period NOT available in dataset)\n\n")
  
  if (ANALYTIC_TIER$permit_rcs) {
    models$C <- cph(Surv(Time_survival, Outcome_OS) ~ 
                      N_stage + rcs(LN_total, 3) + AG + Sex_f + 
                      T_stage_group + Grade_f + LVI_f,
                    data = df, x = TRUE, y = TRUE, surv = TRUE)
  } else {
    models$C <- cph(Surv(Time_survival, Outcome_OS) ~ 
                      N_stage + LN_total + AG + Sex_f + 
                      T_stage_group + Grade_f + LVI_f,
                    data = df, x = TRUE, y = TRUE, surv = TRUE)
  }
  print(models$C)
  cat("\nAnova:\n")
  print(anova(models$C))
  
  model_summary <- rbind(model_summary, data.frame(
    Model = "C", Layer = "3 (Full MVA)",
    Description = "N-stage + all covariates (PRIMARY)",
    df_used = models$C$stats["d.f."],
    AIC = AIC(models$C), BIC = BIC(models$C),
    LogLik = models$C$loglik[2]
  ))
  
  # --- Model D: LNR + TLE + all covariates (Primary LNR model) ---
  cat("\n━━━ Model D: LNR Full MVA (PRIMARY) ━━━\n")
  cat("  Covariates: LNR + TLE + Age + Sex + T-stage + Grade + LVI\n\n")
  
  if (ANALYTIC_TIER$permit_rcs) {
    models$D <- cph(Surv(Time_survival, Outcome_OS) ~ 
                      rcs(LNR_prop, 3) + rcs(LN_total, 3) + AG + Sex_f + 
                      T_stage_group + Grade_f + LVI_f,
                    data = df, x = TRUE, y = TRUE, surv = TRUE)
  } else {
    models$D <- cph(Surv(Time_survival, Outcome_OS) ~ 
                      LNR_prop + LN_total + AG + Sex_f + 
                      T_stage_group + Grade_f + LVI_f,
                    data = df, x = TRUE, y = TRUE, surv = TRUE)
  }
  print(models$D)
  cat("\nAnova:\n")
  print(anova(models$D))
  
  model_summary <- rbind(model_summary, data.frame(
    Model = "D", Layer = "3 (Full MVA)",
    Description = "LNR + all covariates (PRIMARY)",
    df_used = models$D$stats["d.f."],
    AIC = AIC(models$D), BIC = BIC(models$D),
    LogLik = models$D$loglik[2]
  ))
  
} else {
  cat("\n  Layer 3 SKIPPED — MVA not permitted at Tier", 
      ANALYTIC_TIER$tier_code, "\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 4. INTERACTION MODEL (Only if Tier >= 4, i.e. >= 200 events)
# ═══════════════════════════════════════════════════════════════════════════════

if (ANALYTIC_TIER$permit_interaction) {
  section_header("INTERACTION MODEL: LNR × TLE ADEQUACY (H4)")
  
  if (ANALYTIC_TIER$permit_rcs) {
    models$D_interaction <- cph(
      Surv(Time_survival, Outcome_OS) ~ 
        rcs(LNR_prop, 3) * TLE_adequate_binary + rcs(LN_total, 3) + 
        AG + Sex_f + T_stage_group + Grade_f + LVI_f,
      data = df, x = TRUE, y = TRUE, surv = TRUE
    )
  } else {
    models$D_interaction <- cph(
      Surv(Time_survival, Outcome_OS) ~ 
        LNR_prop * TLE_adequate_binary + LN_total + 
        AG + Sex_f + T_stage_group + Grade_f + LVI_f,
      data = df, x = TRUE, y = TRUE, surv = TRUE
    )
  }
  print(models$D_interaction)
  cat("\nAnova (interaction test):\n")
  print(anova(models$D_interaction))
  
} else {
  cat("\n  Interaction model NOT permitted at Tier", 
      ANALYTIC_TIER$tier_code, "(requires >= 200 events)\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 5. PROPORTIONAL HAZARDS ASSUMPTION TESTING
# ═══════════════════════════════════════════════════════════════════════════════

section_header("PROPORTIONAL HAZARDS ASSUMPTION TESTING")

# Test PH for all built models using survival::cox.zph on coxph equivalent
# cph objects from rms can be tested similarly
for (mname in names(models)) {
  cat(sprintf("\n━━━ PH test for Model %s ━━━\n", mname))
  
  # Convert to coxph if needed for cox.zph
  tryCatch({
    ph_test <- cox.zph(models[[mname]], transform = "log")
    print(ph_test)
    
    # Check for violations
    violations <- ph_test$table[, "p"] < PH_VIOLATION_THRESHOLD
    if (any(violations[-length(violations)])) {  # Exclude GLOBAL
      violated_vars <- rownames(ph_test$table)[which(violations[-length(violations)])]
      cat("  *** PH VIOLATION detected for:", 
          paste(violated_vars, collapse = ", "), "\n")
      cat("  --> Per protocol: consider stratification or time-interaction term.\n")
    } else {
      cat("  No PH violations detected (all p >", PH_VIOLATION_THRESHOLD, ").\n")
    }
    
    # Save Schoenfeld residual plots for primary models (C and D)
    if (mname %in% c("C", "D")) {
      pdf(file.path(DIR_FIGURES, paste0("Fig_PH_test_Model_", mname, ".pdf")),
          width = FIG_WIDTH, height = FIG_HEIGHT)
      plot(ph_test, main = paste("Schoenfeld Residuals — Model", mname))
      dev.off()
      cat("  Saved: Fig_PH_test_Model_", mname, ".pdf\n")
    }
  }, error = function(e) {
    cat("  PH test error:", e$message, "\n")
  })
}

# ═══════════════════════════════════════════════════════════════════════════════
# 6. MODEL COMPARISON SUMMARY TABLE
# ═══════════════════════════════════════════════════════════════════════════════

section_header("MODEL COMPARISON SUMMARY")

cat("Model summary across all layers:\n")
print(model_summary)
save_table(model_summary, "model_comparison_summary")

# ═══════════════════════════════════════════════════════════════════════════════
# 7. ASSOCIATION PLOT — LOG-HR VS LNR (FROM RCS MODEL)
# ═══════════════════════════════════════════════════════════════════════════════

section_header("ASSOCIATION PLOTS")

# From the primary Model D (if available)
if (!is.null(models$D)) {
  pred_lnr <- Predict(models$D, LNR_prop, fun = exp)
  
  p_assoc <- ggplot(pred_lnr, aes(x = LNR_prop, y = yhat)) +
    geom_line(colour = "#B2182B", linewidth = 1.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), 
                fill = "#B2182B", alpha = 0.12) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = c(LNR_CUT_LOW_INTER, LNR_CUT_INTER_HIGH),
               linetype = "dotted", colour = c("#4DAF4A", "#FF7F00"),
               linewidth = 0.6) +
    labs(title = "LNR–Mortality Association (Full MVA, Model D)",
         subtitle = "Hazard ratio relative to reference | Adjusted for all covariates",
         x = "Lymph Node Ratio (proportion)",
         y = "Hazard Ratio (HR)") +
    theme_publication()
  save_figure(p_assoc, "Fig_LNR_HR_association_fullMVA")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 8. FOREST PLOT — FULL MVA HAZARD RATIOS
# ═══════════════════════════════════════════════════════════════════════════════

if (!is.null(models$D)) {
  section_header("FOREST PLOT: MODEL D (LNR FULL MVA)")
  
  # Extract coefficients from coxph-equivalent for forest plot
  cox_D <- coxph(Surv(Time_survival, Outcome_OS) ~ 
                   LNR_prop + LN_total + AG + Sex_f + 
                   T_stage_group + Grade_f + LVI_f,
                 data = df)
  
  forest_data <- broom::tidy(cox_D, exponentiate = TRUE, conf.int = TRUE)
  forest_data$term <- gsub("Sex_f", "Sex: ", forest_data$term)
  forest_data$term <- gsub("T_stage_group", "T-stage: ", forest_data$term)
  forest_data$term <- gsub("Grade_f", "Grade: ", forest_data$term)
  forest_data$term <- gsub("LVI_f", "LVI: ", forest_data$term)
  
  p_forest <- ggplot(forest_data, aes(x = estimate, y = reorder(term, estimate))) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_point(size = 3, colour = "#B2182B") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                   height = 0.2, colour = "#B2182B") +
    labs(title = "Forest Plot: Model D (LNR Full MVA)",
         subtitle = "Hazard ratios with 95% CI for Overall Survival",
         x = "Hazard Ratio (log scale)",
         y = NULL) +
    scale_x_log10() +
    theme_publication() +
    theme(axis.text.y = element_text(size = 10))
  save_figure(p_forest, "Fig_forest_plot_ModelD", height = 5)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 9. SAVE ALL MODEL OBJECTS
# ═══════════════════════════════════════════════════════════════════════════════

save(models, model_summary, dd,
     file = file.path(DIR_MODELS, "phase3_cox_models.RData"))

cat("\n", strrep("=", 78), "\n")
cat("  SCRIPT 04 COMPLETE\n")
cat("  Models built:", paste(names(models), collapse = ", "), "\n")
cat(strrep("=", 78), "\n")
