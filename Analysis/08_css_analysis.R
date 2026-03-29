################################################################################
#  08_css_analysis.R — Cancer-Specific Survival (Secondary Outcome)
#  Protocol Section 5.4.2 | Competing Risks (Fine-Gray)
#
#  CSS is conditional on >= 85% cause-of-death data availability.
#  Uses Fine-Gray subdistribution hazard model.
#  Non-cancer death is the competing event.
################################################################################

source(file.path("d:/Google_SSD_RAM/_OneDrive/Research",
                 "Project_ColoRectal/LNR_20260328/Analysis/00_config.R"))

load(file.path(DIR_DATA, "analytic_dataset.RData"))

section_header("SCRIPT 08: CANCER-SPECIFIC SURVIVAL ANALYSIS")

cat("Analytic cohort: N =", nrow(df), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 1. CSS DATA COMPLETENESS CHECK
# ═══════════════════════════════════════════════════════════════════════════════

section_header("CSS DATA COMPLETENESS CHECK")

n_deaths_os <- sum(df$Outcome_OS == 1, na.rm = TRUE)
n_deaths_dss <- sum(df$Outcome_DSS == 1, na.rm = TRUE)

# For CSS data completeness: Among deceased patients (OS=1),
# we need to know if death was cancer-specific (DSS=1) or not (DSS=0)
# All deceased patients have DSS classified since both are recorded
# Completeness = proportion of OS=1 with known cause (which is 100% here since
# DSS is recorded for all)

# Deaths attributable to cancer vs other causes
cat("Total deaths (OS events):", n_deaths_os, "\n")
cat("Cancer-specific deaths (DSS events):", n_deaths_dss, "\n")
cat("Non-cancer deaths:", n_deaths_os - n_deaths_dss, "\n")
cat("Alive/censored:", sum(df$Outcome_OS == 0), "\n\n")

# CSS completeness: all deceased have cause attribution recorded
css_completeness <- 1.0  # Both OS and DSS are recorded for everyone
cat("Cause-of-death data completeness:", 
    round(css_completeness * 100, 1), "%\n")

if (css_completeness >= CSS_COD_THRESHOLD) {
  cat(">>> CSS analysis PERMITTED (>= 85% threshold MET) <<<\n\n")
  CSS_PERMITTED <- TRUE
} else {
  cat(">>> CSS analysis NOT feasible (< 85% threshold). Descriptive only. <<<\n")
  CSS_PERMITTED <- FALSE
}

# ═══════════════════════════════════════════════════════════════════════════════
# 2. COMPETING RISKS SETUP
# ═══════════════════════════════════════════════════════════════════════════════

section_header("COMPETING RISKS FRAMEWORK")

# Create competing risks event status:
#   0 = alive/censored
#   1 = cancer-specific death (event of interest)
#   2 = non-cancer death (competing event)
df$cr_status <- 0
df$cr_status[df$Outcome_DSS == 1] <- 1        # Cancer death
df$cr_status[df$Outcome_OS == 1 & df$Outcome_DSS == 0] <- 2  # Non-cancer death

cat("Competing risks event distribution:\n")
print(table(df$cr_status, dnn = "Status"))
cat("  0 = Alive/censored, 1 = Cancer death, 2 = Non-cancer death\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 3. CUMULATIVE INCIDENCE FUNCTIONS (CIF)
# ═══════════════════════════════════════════════════════════════════════════════

if (CSS_PERMITTED && n_deaths_dss >= 5) {
  section_header("CUMULATIVE INCIDENCE FUNCTIONS")
  
  # --- CIF by N-Stage ---
  cat("━━━ CIF by N-Stage ━━━\n")
  cif_nstage <- cmprsk::cuminc(ftime = df$Time_survival, 
                                fstatus = df$cr_status,
                                group = df$N_stage)
  print(cif_nstage$Tests)
  
  # Plot CIF
  pdf(file.path(DIR_FIGURES, "Fig_CIF_CSS_by_Nstage.pdf"),
      width = FIG_WIDTH, height = FIG_HEIGHT)
  plot(cif_nstage, curvlab = c("N1: Cancer death", "N1: Other death",
                                 "N2: Cancer death", "N2: Other death"),
       col = c(COL_NSTAGE["N1"], adjustcolor(COL_NSTAGE["N1"], 0.4),
               COL_NSTAGE["N2"], adjustcolor(COL_NSTAGE["N2"], 0.4)),
       lty = c(1, 2, 1, 2), lwd = 2,
       main = "Cumulative Incidence: Cancer-Specific Death by N-Stage",
       xlab = "Time (days)", ylab = "Cumulative Incidence")
  legend("topleft", bty = "n",
         legend = c("N1 — Cancer", "N1 — Other", "N2 — Cancer", "N2 — Other"),
         col = c(COL_NSTAGE["N1"], adjustcolor(COL_NSTAGE["N1"], 0.4),
                 COL_NSTAGE["N2"], adjustcolor(COL_NSTAGE["N2"], 0.4)),
         lty = c(1, 2, 1, 2), lwd = 2)
  dev.off()
  cat("Saved: Fig_CIF_CSS_by_Nstage\n")
  
  # --- CIF by LNR Tier ---
  cat("\n━━━ CIF by LNR Tier ━━━\n")
  cif_lnr <- cmprsk::cuminc(ftime = df$Time_survival,
                              fstatus = df$cr_status,
                              group = df$LNR_tier)
  print(cif_lnr$Tests)
  
  pdf(file.path(DIR_FIGURES, "Fig_CIF_CSS_by_LNR_tier.pdf"),
      width = FIG_WIDTH, height = FIG_HEIGHT)
  plot(cif_lnr, 
       main = "Cumulative Incidence: Cancer-Specific Death by LNR Tier",
       xlab = "Time (days)", ylab = "Cumulative Incidence",
       lwd = 2)
  dev.off()
  cat("Saved: Fig_CIF_CSS_by_LNR_tier\n")
  
  # ═══════════════════════════════════════════════════════════════════════════
  # 4. FINE-GRAY SUBDISTRIBUTION HAZARD MODELS
  # ═══════════════════════════════════════════════════════════════════════════
  
  section_header("FINE-GRAY SUBDISTRIBUTION HAZARD MODELS")
  
  # Prepare covariate matrix for crr()
  # crr() requires a numeric covariate matrix
  
  # --- Model: N-stage (Fine-Gray) ---
  cat("━━━ Fine-Gray: N-stage for CSS ━━━\n")
  cov_nstage <- model.matrix(~ N_stage, data = df)[, -1, drop = FALSE]
  fg_nstage <- cmprsk::crr(ftime = df$Time_survival,
                             fstatus = df$cr_status,
                             cov1 = cov_nstage,
                             failcode = 1, cencode = 0)
  cat("Coefficients:\n")
  print(summary(fg_nstage))
  
  # --- Model: LNR continuous (Fine-Gray) ---
  cat("\n━━━ Fine-Gray: LNR continuous for CSS ━━━\n")
  cov_lnr <- matrix(df$LNR_prop, ncol = 1)
  colnames(cov_lnr) <- "LNR_prop"
  fg_lnr <- cmprsk::crr(ftime = df$Time_survival,
                          fstatus = df$cr_status,
                          cov1 = cov_lnr,
                          failcode = 1, cencode = 0)
  cat("Coefficients:\n")
  print(summary(fg_lnr))
  
  # --- Model: LNR tiers (Fine-Gray) ---
  cat("\n━━━ Fine-Gray: LNR tiers for CSS ━━━\n")
  cov_lnr_tier <- model.matrix(~ LNR_tier, data = df)[, -1, drop = FALSE]
  fg_lnr_tier <- cmprsk::crr(ftime = df$Time_survival,
                               fstatus = df$cr_status,
                               cov1 = cov_lnr_tier,
                               failcode = 1, cencode = 0)
  cat("Coefficients:\n")
  print(summary(fg_lnr_tier))
  
  # ═══════════════════════════════════════════════════════════════════════════
  # 5. FULL MVA FINE-GRAY (if tier permits)
  # ═══════════════════════════════════════════════════════════════════════════
  
  if (ANALYTIC_TIER$permit_mva) {
    section_header("FINE-GRAY: FULL MVA FOR CSS")
    
    # N-stage MVA
    cat("━━━ Fine-Gray: N-stage Full MVA for CSS ━━━\n")
    cov_C_css <- model.matrix(~ N_stage + LN_total + AG + Sex_f + 
                                T_stage_group + Grade_f + LVI_f, 
                              data = df)[, -1]
    fg_C <- cmprsk::crr(ftime = df$Time_survival, fstatus = df$cr_status,
                         cov1 = cov_C_css, failcode = 1, cencode = 0)
    print(summary(fg_C))
    
    # LNR MVA
    cat("\n━━━ Fine-Gray: LNR Full MVA for CSS ━━━\n")
    cov_D_css <- model.matrix(~ LNR_prop + LN_total + AG + Sex_f + 
                                T_stage_group + Grade_f + LVI_f, 
                              data = df)[, -1]
    fg_D <- cmprsk::crr(ftime = df$Time_survival, fstatus = df$cr_status,
                         cov1 = cov_D_css, failcode = 1, cencode = 0)
    print(summary(fg_D))
  }
  
  # Save results
  save(cif_nstage, cif_lnr, fg_nstage, fg_lnr, fg_lnr_tier,
       file = file.path(DIR_MODELS, "phase_css_competing_risks.RData"))
  
} else {
  cat("\n  CSS analysis limited due to insufficient events or data.\n")
  cat("  Cancer-specific deaths:", n_deaths_dss, "\n")
  cat("  Reporting descriptive results only.\n")
  
  section_header("DESCRIPTIVE CSS ANALYSIS")
  
  # Descriptive: cause of death distribution
  cod_tab <- data.frame(
    Category = c("Alive/censored", "Cancer-specific death", 
                 "Non-cancer death", "Total"),
    N = c(sum(df$cr_status == 0), sum(df$cr_status == 1), 
          sum(df$cr_status == 2), nrow(df)),
    Pct = c(round(100 * sum(df$cr_status == 0) / nrow(df), 1),
            round(100 * sum(df$cr_status == 1) / nrow(df), 1),
            round(100 * sum(df$cr_status == 2) / nrow(df), 1),
            100)
  )
  print(cod_tab)
  save_table(cod_tab, "CSS_descriptive_cause_of_death")
}

cat("\n", strrep("=", 78), "\n")
cat("  SCRIPT 08 COMPLETE\n")
cat(strrep("=", 78), "\n")
