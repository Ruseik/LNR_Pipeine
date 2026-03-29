################################################################################
#  05_performance_metrics.R — Phase 3: Performance Metrics
#  Protocol Section 7.4 | SO2
#
#  Primary co-endpoints: Uno C-statistic + Decision-Curve Analysis
#  Primary secondary: Calibration + AIC/BIC
#  Exploratory: NRI, IDI (appendix only)
#  Internal validation: Bootstrap optimism-corrected C-statistic
################################################################################

source(file.path("d:/Google_SSD_RAM/_OneDrive/Research",
                 "Project_ColoRectal/LNR_20260328/Analysis/00_config.R"))

load(file.path(DIR_DATA, "analytic_dataset.RData"))
load(file.path(DIR_MODELS, "phase3_cox_models.RData"))

section_header("SCRIPT 05: PERFORMANCE METRICS (Phase 3)")

cat("Analytic tier:", ANALYTIC_TIER$label, "\n\n")

# Re-establish rms datadist
dd <- datadist(df)
options(datadist = "dd")

# ═══════════════════════════════════════════════════════════════════════════════
# 1. UNO'S C-STATISTIC FOR ALL MODELS
# ═══════════════════════════════════════════════════════════════════════════════

section_header("UNO'S C-STATISTIC")

# Truncation time for Uno's C
tau_uno <- quantile(df$Time_survival[df$Outcome_OS == 1], 
                    UNO_TRUNCATION_QUANTILE, na.rm = TRUE)
cat("Uno C-stat truncation (75th percentile of event times):", 
    round(tau_uno, 0), "days\n\n")

# Function to compute Uno C with bootstrap CI
compute_uno_c <- function(model_name, model_obj, data) {
  # Get linear predictor
  lp <- tryCatch(predict(model_obj, newdata = data, type = "lp"),
                 error = function(e) predict(model_obj, type = "lp"))
  
  # Uno's concordance via survival::concordance
  conc <- concordance(Surv(Time_survival, Outcome_OS) ~ lp, 
                      data = data, timewt = "n/G2")
  
  c_stat <- conc$concordance
  c_se <- sqrt(conc$var)
  c_lower <- c_stat - 1.96 * c_se
  c_upper <- c_stat + 1.96 * c_se
  
  data.frame(
    Model = model_name,
    C_statistic = round(c_stat, 4),
    SE = round(c_se, 4),
    CI_lower = round(c_lower, 4),
    CI_upper = round(c_upper, 4),
    stringsAsFactors = FALSE
  )
}

# Compute C for all models
c_results <- data.frame()
for (mname in names(models)) {
  tryCatch({
    c_row <- compute_uno_c(mname, models[[mname]], df)
    c_results <- rbind(c_results, c_row)
    cat(sprintf("  Model %s: C = %.4f (95%% CI: %.4f - %.4f)\n",
                c_row$Model, c_row$C_statistic, c_row$CI_lower, c_row$CI_upper))
  }, error = function(e) {
    cat(sprintf("  Model %s: C-stat computation failed: %s\n", mname, e$message))
  })
}

cat("\nC-statistic summary table:\n")
print(c_results)
save_table(c_results, "C_statistic_all_models")

# ═══════════════════════════════════════════════════════════════════════════════
# 2. PRIMARY COMPARISON: ΔC BETWEEN MODEL D AND MODEL C
# ═══════════════════════════════════════════════════════════════════════════════

section_header("PRIMARY ΔC COMPARISON (MODEL D vs MODEL C)")

if (!is.null(models$C) && !is.null(models$D)) {
  
  c_model_C <- c_results$C_statistic[c_results$Model == "C"]
  c_model_D <- c_results$C_statistic[c_results$Model == "D"]
  delta_c <- c_model_D - c_model_C
  
  cat("Model C (N-stage Full MVA) C-statistic:", c_model_C, "\n")
  cat("Model D (LNR Full MVA) C-statistic:", c_model_D, "\n")
  cat("ΔC (Model D - Model C):", round(delta_c, 4), "\n\n")
  
  # Bootstrap ΔC with BCa CI
  cat("Computing bootstrap 95% CI for ΔC (", N_BOOTSTRAP, " resamples)...\n")
  
  boot_delta_c <- function(data, indices) {
    boot_data <- data[indices, ]
    
    # Refit both models on bootstrap sample
    tryCatch({
      if (ANALYTIC_TIER$permit_rcs) {
        dd_boot <- datadist(boot_data)
        options(datadist = "dd_boot")
        
        fit_C <- cph(Surv(Time_survival, Outcome_OS) ~ 
                       N_stage + rcs(LN_total, 3) + AG + Sex_f + 
                       T_stage_group + Grade_f + LVI_f,
                     data = boot_data, x = TRUE, y = TRUE)
        fit_D <- cph(Surv(Time_survival, Outcome_OS) ~ 
                       rcs(LNR_prop, 3) + rcs(LN_total, 3) + AG + Sex_f + 
                       T_stage_group + Grade_f + LVI_f,
                     data = boot_data, x = TRUE, y = TRUE)
      } else {
        fit_C <- cph(Surv(Time_survival, Outcome_OS) ~ 
                       N_stage + LN_total + AG + Sex_f + 
                       T_stage_group + Grade_f + LVI_f,
                     data = boot_data, x = TRUE, y = TRUE)
        fit_D <- cph(Surv(Time_survival, Outcome_OS) ~ 
                       LNR_prop + LN_total + AG + Sex_f + 
                       T_stage_group + Grade_f + LVI_f,
                     data = boot_data, x = TRUE, y = TRUE)
      }
      
      lp_C <- predict(fit_C, type = "lp")
      lp_D <- predict(fit_D, type = "lp")
      
      conc_C <- concordance(Surv(Time_survival, Outcome_OS) ~ lp_C, 
                            data = boot_data, timewt = "n/G2")
      conc_D <- concordance(Surv(Time_survival, Outcome_OS) ~ lp_D,
                            data = boot_data, timewt = "n/G2")
      
      return(conc_D$concordance - conc_C$concordance)
    }, error = function(e) {
      return(NA)
    })
  }
  
  set.seed(SEED)
  boot_result <- boot(data = df, statistic = boot_delta_c, R = N_BOOTSTRAP)
  
  # BCa confidence interval
  boot_ci <- tryCatch(
    boot.ci(boot_result, type = "bca", conf = 1 - ALPHA),
    error = function(e) {
      cat("  BCa CI failed, using percentile method.\n")
      boot.ci(boot_result, type = "perc", conf = 1 - ALPHA)
    }
  )
  
  if (!is.null(boot_ci$bca)) {
    delta_c_ci <- boot_ci$bca[4:5]
  } else if (!is.null(boot_ci$percent)) {
    delta_c_ci <- boot_ci$percent[4:5]
  } else {
    delta_c_ci <- quantile(boot_result$t, c(0.025, 0.975), na.rm = TRUE)
  }
  
  cat("\n═══ PRIMARY ESTIMAND RESULT ═══\n")
  cat("ΔC (Model D - Model C):", round(delta_c, 4), "\n")
  cat("95% Bootstrap CI: [", round(delta_c_ci[1], 4), ",", 
      round(delta_c_ci[2], 4), "]\n")
  cat("Pre-specified threshold (ΔC ≥ 0.05):",
      ifelse(delta_c >= DELTA_C_THRESHOLD, "MET", "NOT MET"), "\n")
  cat("CI excludes zero:", 
      ifelse(delta_c_ci[1] > 0, "YES (Model D superior)", 
             ifelse(delta_c_ci[2] < 0, "YES (Model C superior)", "NO")), "\n")
  
  # Save primary result
  primary_result <- data.frame(
    Comparison = "Model D vs Model C (Layer 3)",
    Delta_C = round(delta_c, 4),
    CI_lower = round(delta_c_ci[1], 4),
    CI_upper = round(delta_c_ci[2], 4),
    Threshold_met = delta_c >= DELTA_C_THRESHOLD,
    Conclusion = ifelse(delta_c >= DELTA_C_THRESHOLD & delta_c_ci[1] > 0,
                        "LNR significantly superior",
                        ifelse(delta_c_ci[2] < 0,
                               "N-stage superior",
                               "No significant difference"))
  )
  save_table(primary_result, "PRIMARY_delta_C_result")
  
  # Also compute ΔC for Layer 1 and Layer 2 to show the trend
  if (!is.null(models$A) && !is.null(models$B)) {
    c_A <- c_results$C_statistic[c_results$Model == "A"]
    c_B <- c_results$C_statistic[c_results$Model == "B"]
    cat("\nLayer 1 ΔC (B - A):", round(c_B - c_A, 4), "\n")
  }
  if (!is.null(models$C2) && !is.null(models$D2)) {
    c_C2 <- c_results$C_statistic[c_results$Model == "C2"]
    c_D2 <- c_results$C_statistic[c_results$Model == "D2"]
    cat("Layer 2 ΔC (D2 - C2):", round(c_D2 - c_C2, 4), "\n")
  }

} else {
  cat("  Models C and D not available (tier too low for MVA).\n")
  cat("  Using Layer 1 comparison: Model B vs Model A.\n")
  
  if (!is.null(models$A) && !is.null(models$B)) {
    c_A <- c_results$C_statistic[c_results$Model == "A"]
    c_B <- c_results$C_statistic[c_results$Model == "B"]
    delta_c <- c_B - c_A
    cat("ΔC (Model B - Model A):", round(delta_c, 4), "\n")
    cat("NOTE: This is a hypothesis-generating comparison (inadequate tier).\n")
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# 3. DECISION-CURVE ANALYSIS (Co-Primary)
# ═══════════════════════════════════════════════════════════════════════════════

if (ANALYTIC_TIER$permit_dca) {
  section_header("DECISION-CURVE ANALYSIS (Co-Primary)")
  
  # DCA requires predicted probabilities of the event
  # Use 5-year survival prediction as the basis
  
  # Predicted probabilities from coxph models
  # Using 5-year risk = 1 - survival probability at 5 years
  
  tryCatch({
    # Refit with coxph for easier prediction
    if (ANALYTIC_TIER$permit_rcs) {
      cox_C_dca <- coxph(Surv(Time_survival, Outcome_OS) ~ 
                           N_stage + rcs(LN_total, 3) + AG + Sex_f + 
                           T_stage_group + Grade_f + LVI_f,
                         data = df)
      cox_D_dca <- coxph(Surv(Time_survival, Outcome_OS) ~ 
                           rcs(LNR_prop, 3) + rcs(LN_total, 3) + AG + Sex_f + 
                           T_stage_group + Grade_f + LVI_f,
                         data = df)
    } else {
      cox_C_dca <- coxph(Surv(Time_survival, Outcome_OS) ~ 
                           N_stage + LN_total + AG + Sex_f + 
                           T_stage_group + Grade_f + LVI_f,
                         data = df)
      cox_D_dca <- coxph(Surv(Time_survival, Outcome_OS) ~ 
                           LNR_prop + LN_total + AG + Sex_f + 
                           T_stage_group + Grade_f + LVI_f,
                         data = df)
    }
    
    # Baseline survival at 5 years (1826 days)
    t_star <- RMST_TAU_DAYS
    baseline_C <- basehaz(cox_C_dca, centered = FALSE)
    baseline_D <- basehaz(cox_D_dca, centered = FALSE)
    
    # 5-year predicted risk for each patient
    H0_C <- approx(baseline_C$time, baseline_C$hazard, xout = t_star, rule = 2)$y
    H0_D <- approx(baseline_D$time, baseline_D$hazard, xout = t_star, rule = 2)$y
    
    lp_C <- predict(cox_C_dca, type = "lp")
    lp_D <- predict(cox_D_dca, type = "lp")
    
    df$risk_C <- 1 - exp(-H0_C * exp(lp_C))
    df$risk_D <- 1 - exp(-H0_D * exp(lp_D))
    
    # Create binary outcome for 5-year event
    df$event_5yr <- ifelse(df$Outcome_OS == 1 & df$Time_survival <= t_star, 1, 0)
    
    # Run DCA
    dca_result <- dca(
      Surv(Time_survival, Outcome_OS) ~ risk_C + risk_D,
      data = df,
      time = t_star,
      thresholds = DCA_THRESHOLDS
    )
    
    # DCA plot
    p_dca <- plot(dca_result) +
      labs(title = "Decision Curve Analysis: Model C vs Model D",
           subtitle = paste0("5-year overall survival | Threshold range: ",
                             DCA_THRESHOLD_MIN * 100, "%-", 
                             DCA_THRESHOLD_MAX * 100, "%")) +
      theme_publication()
    save_figure(p_dca, "Fig_DCA_ModelC_vs_ModelD")
    
    cat("DCA plot saved. Visual interpretation:\n")
    cat("  - Model with higher net benefit across threshold range is preferred.\n")
    cat("  - If curves are close, neither model has clear clinical advantage.\n")
    
  }, error = function(e) {
    cat("DCA computation encountered an error:", e$message, "\n")
    cat("This may occur with small sample sizes or sparse events.\n")
  })
  
} else {
  cat("\n  DCA not permitted at Tier", ANALYTIC_TIER$tier_code, "\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 4. CALIBRATION (Primary Secondary Endpoint)
# ═══════════════════════════════════════════════════════════════════════════════

if (ANALYTIC_TIER$permit_mva && !is.null(models$D)) {
  section_header("CALIBRATION ASSESSMENT")
  
  for (mname in c("C", "D")) {
    if (!is.null(models[[mname]])) {
      cat(sprintf("\n━━━ Calibration: Model %s ━━━\n", mname))
      
      tryCatch({
        cal <- calibrate(models[[mname]], B = min(N_BOOTSTRAP, 200), 
                         u = RMST_TAU_DAYS, cmethod = "KM")
        
        pdf(file.path(DIR_FIGURES, paste0("Fig_calibration_Model_", mname, ".pdf")),
            width = FIG_WIDTH, height = FIG_HEIGHT)
        plot(cal, main = paste("Calibration Plot — Model", mname,
                               "(5-year OS)"),
             xlab = "Predicted 5-Year Survival",
             ylab = "Observed 5-Year Survival",
             subtitles = TRUE)
        abline(0, 1, lty = 2, col = "grey50")
        dev.off()
        
        cat("  Calibration plot saved.\n")
      }, error = function(e) {
        cat("  Calibration failed:", e$message, "\n")
        cat("  This may be due to insufficient events or extreme predictions.\n")
      })
    }
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# 5. BOOTSTRAP OPTIMISM-CORRECTED C-STATISTIC (Internal Validation)
# ═══════════════════════════════════════════════════════════════════════════════

if (ANALYTIC_TIER$permit_mva) {
  section_header("BOOTSTRAP INTERNAL VALIDATION")
  
  for (mname in c("C", "D")) {
    if (!is.null(models[[mname]])) {
      cat(sprintf("\n━━━ Bootstrap validation: Model %s ━━━\n", mname))
      
      tryCatch({
        val <- validate(models[[mname]], B = N_BOOTSTRAP, method = "boot")
        print(val)
        
        # Extract optimism-corrected Dxy → C
        dxy_corrected <- val["Dxy", "index.corrected"]
        c_corrected <- 0.5 * (dxy_corrected + 1)
        cat(sprintf("  Optimism-corrected C-statistic: %.4f\n", c_corrected))
        cat(sprintf("  Optimism: %.4f\n", val["Dxy", "optimism"] / 2))
        
      }, error = function(e) {
        cat("  Bootstrap validation failed:", e$message, "\n")
      })
    }
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# 6. AIC / BIC COMPARISON
# ═══════════════════════════════════════════════════════════════════════════════

section_header("AIC / BIC COMPARISON")

aic_bic <- model_summary[, c("Model", "Layer", "AIC", "BIC")]
cat("AIC / BIC comparison across models:\n")
print(aic_bic)

# Identify best model by AIC
best_aic <- aic_bic$Model[which.min(aic_bic$AIC)]
best_bic <- aic_bic$Model[which.min(aic_bic$BIC)]
cat("\nBest model by AIC:", best_aic, "\n")
cat("Best model by BIC:", best_bic, "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 7. EXPLORATORY: NRI AND IDI (Appendix Only)
# ═══════════════════════════════════════════════════════════════════════════════

if (ANALYTIC_TIER$permit_mva && !is.null(models$C) && !is.null(models$D)) {
  section_header("EXPLORATORY: NRI AND IDI (APPENDIX ONLY)")
  cat("NOTE: These metrics are reported in supplementary appendix only.\n")
  cat("They are NOT used for primary conclusions per protocol.\n\n")
  
  tryCatch({
    # Predicted probabilities (linear predictors)
    lp_model_C <- predict(models$C, type = "lp")
    lp_model_D <- predict(models$D, type = "lp")
    
    # Continuous NRI using survIDINRI or nricens package
    # Simple implementation: direction of reclassification
    improved <- (lp_model_D > lp_model_C & df$Outcome_OS == 1) |
                (lp_model_D < lp_model_C & df$Outcome_OS == 0)
    worsened <- (lp_model_D < lp_model_C & df$Outcome_OS == 1) |
                (lp_model_D > lp_model_C & df$Outcome_OS == 0)
    
    nri_events <- mean(lp_model_D[df$Outcome_OS == 1] > 
                        lp_model_C[df$Outcome_OS == 1]) -
                  mean(lp_model_D[df$Outcome_OS == 1] < 
                        lp_model_C[df$Outcome_OS == 1])
    nri_nonevents <- mean(lp_model_D[df$Outcome_OS == 0] < 
                           lp_model_C[df$Outcome_OS == 0]) -
                     mean(lp_model_D[df$Outcome_OS == 0] > 
                           lp_model_C[df$Outcome_OS == 0])
    nri_continuous <- nri_events + nri_nonevents
    
    cat("Continuous NRI (Pencina 2011):\n")
    cat("  NRI_events:", round(nri_events, 4), "\n")
    cat("  NRI_nonevents:", round(nri_nonevents, 4), "\n")
    cat("  NRI_total:", round(nri_continuous, 4), "\n\n")
    
    # IDI
    idi <- mean(lp_model_D[df$Outcome_OS == 1]) - 
           mean(lp_model_D[df$Outcome_OS == 0]) -
           (mean(lp_model_C[df$Outcome_OS == 1]) - 
            mean(lp_model_C[df$Outcome_OS == 0]))
    cat("IDI:", round(idi, 4), "\n")
    
    nri_idi_result <- data.frame(
      Metric = c("NRI_events", "NRI_nonevents", "NRI_total", "IDI"),
      Value = round(c(nri_events, nri_nonevents, nri_continuous, idi), 4),
      Note = rep("Exploratory — appendix only", 4)
    )
    save_table(nri_idi_result, "appendix_NRI_IDI")
    
  }, error = function(e) {
    cat("NRI/IDI computation failed:", e$message, "\n")
  })
}

# ═══════════════════════════════════════════════════════════════════════════════
# 8. C-STATISTIC COMPARISON PLOT
# ═══════════════════════════════════════════════════════════════════════════════

section_header("C-STATISTIC COMPARISON PLOT")

if (nrow(c_results) > 0) {
  c_results$Model_label <- paste0("Model ", c_results$Model)
  c_results$Model_label <- factor(c_results$Model_label, 
                                   levels = rev(c_results$Model_label))
  
  # Assign layers for colouring
  c_results$Layer <- ifelse(c_results$Model %in% c("A", "B"), "Layer 1",
                     ifelse(c_results$Model %in% c("C2", "D2"), "Layer 2",
                            "Layer 3"))
  
  p_cstat <- ggplot(c_results, aes(x = C_statistic, y = Model_label, 
                                    colour = Layer)) +
    geom_point(size = 4) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.3) +
    geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey60") +
    scale_colour_manual(values = c("Layer 1" = "#4393C3", 
                                    "Layer 2" = "#2166AC",
                                    "Layer 3" = "#053061")) +
    labs(title = "C-Statistic Comparison Across Model Layers",
         subtitle = "Uno's C-statistic with 95% CI",
         x = "C-statistic",
         y = NULL,
         colour = "Model Layer") +
    theme_publication() +
    theme(legend.position = "right")
  save_figure(p_cstat, "Fig_C_statistic_comparison")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 9. SAVE RESULTS
# ═══════════════════════════════════════════════════════════════════════════════

save(c_results, model_summary,
     file = file.path(DIR_MODELS, "phase3_performance_metrics.RData"))

cat("\n", strrep("=", 78), "\n")
cat("  SCRIPT 05 COMPLETE\n")
cat(strrep("=", 78), "\n")
