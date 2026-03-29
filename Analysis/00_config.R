################################################################################
#  00_config.R — Master Configuration for LNR vs N-Stage Analysis
#  Protocol: LNR vs TNM N-Stage as a Prognostic Marker in Node-Positive CRC
#  Version 2.0 | Post-Critical-Review Edition
#
#  This script is sourced by ALL subsequent analysis scripts.
#  Modify settings here ONLY — downstream scripts inherit everything.
################################################################################

# ═══════════════════════════════════════════════════════════════════════════════
# 1. PACKAGE MANAGEMENT
# ═══════════════════════════════════════════════════════════════════════════════

# Set up user library path (system library may be read-only)
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

required_packages <- c(
  # -- Core survival analysis
  "survival",       # Cox PH, Surv objects, survfit, cox.zph
  "rms",            # Restricted cubic splines, cph, validate, calibrate
  "survminer",      # Publication-quality KM curves (ggsurvplot)
  
  # -- Competing risks
  "cmprsk",         # Fine-Gray subdistribution hazard models
  "tidycmprsk",     # Tidy interface for cmprsk
  
  # -- Performance metrics
  # Note: Uno's C-statistic is computed via survival::concordance(timewt="n/G2")
  "dcurves",        # Decision-curve analysis
  "pec",            # Prediction error curves (alternative C-stat)
  
  # -- Multiple imputation
  "mice",           # Multiple imputation by chained equations
  
  # -- Tables and reporting
  "tableone",       # Baseline characteristics tables with SMD
  "knitr",          # Table formatting
  "kableExtra",     # Enhanced table output
  "gtsummary",      # Publication-ready summary tables
  

  # -- Graphics
  "ggplot2",        # Core plotting
  "patchwork",      # Multi-panel figure assembly
  "scales",         # Axis formatting
  "RColorBrewer",   # Colour palettes
  "viridis",        # Colourblind-safe palettes
  "ggrepel",        # Non-overlapping text labels
  
  # -- Utilities
  "dplyr",          # Data manipulation
  "tidyr",          # Data reshaping
  "forcats",        # Factor manipulation
  "boot",           # Bootstrap resampling
  "broom"           # Tidy model outputs
)

# Install missing packages
install_if_missing <- function(pkgs) {
  missing <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(missing) > 0) {
    cat("Installing missing packages:", paste(missing, collapse = ", "), "\n")
    install.packages(missing, repos = "https://cloud.r-project.org", 
                     dependencies = TRUE, quiet = TRUE)
  }
}
install_if_missing(required_packages)

# Load all packages
invisible(lapply(required_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

cat("All packages loaded successfully.\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 2. FILE PATHS
# ═══════════════════════════════════════════════════════════════════════════════

# -- Base project directory
DIR_PROJECT <- file.path("d:/Google_SSD_RAM/_OneDrive/Research",
                         "Project_ColoRectal/LNR_20260328")

# -- Analysis scripts directory (this folder)
DIR_ANALYSIS <- file.path(DIR_PROJECT, "Analysis")

# -- Input data
FILE_RAW_DATA <- file.path("d:/Google_SSD_RAM/_OneDrive/Research",
                           "Project_ColoRectal/Analysis_Project_ColoRectal",
                           "CRCALVS_8_1.csv")

# -- Output directories (created if they don't exist)
DIR_OUTPUT    <- file.path(DIR_ANALYSIS, "Outputs")
DIR_TABLES    <- file.path(DIR_OUTPUT, "Tables")
DIR_FIGURES   <- file.path(DIR_OUTPUT, "Figures")
DIR_DATA      <- file.path(DIR_OUTPUT, "Data")
DIR_MODELS    <- file.path(DIR_OUTPUT, "Models")
DIR_APPENDIX  <- file.path(DIR_OUTPUT, "Appendix")

for (d in c(DIR_OUTPUT, DIR_TABLES, DIR_FIGURES, DIR_DATA, DIR_MODELS, DIR_APPENDIX)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 3. REPRODUCIBILITY
# ═══════════════════════════════════════════════════════════════════════════════

SEED <- 2024
set.seed(SEED)

# ═══════════════════════════════════════════════════════════════════════════════
# 4. ANALYTIC CONSTANTS (from Protocol v2.0)
# ═══════════════════════════════════════════════════════════════════════════════

# --- Node-Positive Filtering Strategy ---
# Set to "SN_main" to filter by N-stage (SN_main >= 1)
# Set to "LN_positive" to filter by positive lymph node count (LN_positive >= 1)
# The user can toggle this to compare both approaches
NODE_POSITIVE_FILTER <- "LN_positive"   # <-- ADJUSTABLE: "SN_main" or "LN_positive"

# --- LNR Rosenberg (2008) categorical cutpoints (on proportion scale 0-1) ---
LNR_CUT_LOW_INTER  <- 0.05    # Boundary between Low and Intermediate
LNR_CUT_INTER_HIGH <- 0.20    # Boundary between Intermediate and High
LNR_TIER_LABELS    <- c("Low (<0.05)", "Intermediate (0.05-0.20)", "High (>0.20)")

# --- Alternative LNR cutpoints for sensitivity analyses ---
LNR_ALT_CUT_1 <- c(0.05, 0.40)    # Sensitivity: wider High tier
LNR_ALT_CUT_2 <- 0.20             # Sensitivity: single binary cutpoint

# --- TLE (Total Lymph nodes Examined) adequacy thresholds ---
TLE_ADEQUATE_THRESHOLD <- 12       # AJCC minimum adequacy threshold
TLE_STRINGENT_THRESHOLD <- 8       # Stringent adequacy criterion

# --- RCS (Restricted Cubic Splines) specification ---
RCS_KNOTS_DEFAULT <- 3             # Number of knots for RCS (yields 2 df)
# Knots placed at 10th, 50th, 90th percentiles by default (rms::rcs default)

# --- Uno C-statistic ---
UNO_TRUNCATION_QUANTILE <- 0.75    # Truncation at 75th percentile of event times

# --- Bootstrap ---
N_BOOTSTRAP <- 1000                # Number of bootstrap replicates
BOOTSTRAP_METHOD <- "bca"          # Bias-corrected and accelerated

# --- DCA (Decision Curve Analysis) ---
DCA_THRESHOLD_MIN <- 0.10
DCA_THRESHOLD_MAX <- 0.60
DCA_THRESHOLD_STEP <- 0.01
DCA_THRESHOLDS <- seq(DCA_THRESHOLD_MIN, DCA_THRESHOLD_MAX, by = DCA_THRESHOLD_STEP)

# --- Multiple Imputation ---
MICE_M <- 20                       # Number of imputed datasets
MICE_MISSINGNESS_THRESHOLD <- 0.05 # Apply MICE only if missingness > 5%

# --- Significance and decision thresholds ---
ALPHA <- 0.05                      # Two-sided alpha
DELTA_C_THRESHOLD <- 0.05          # Minimum clinically meaningful ΔC
PH_VIOLATION_THRESHOLD <- 0.10     # Schoenfeld residual p threshold for PH

# --- Landmark analysis ---
LANDMARK_DAYS <- 180               # 6-month landmark (in days)

# --- CSS conditional inclusion ---
CSS_COD_THRESHOLD <- 0.85          # CSS included only if ≥85% cause-of-death data

# --- KM survival time points (years) ---
SURV_TIMEPOINTS_YEARS <- c(1, 3, 5)
SURV_TIMEPOINTS_DAYS  <- SURV_TIMEPOINTS_YEARS * 365.25

# --- RMST truncation ---
RMST_TAU_YEARS <- 5
RMST_TAU_DAYS  <- RMST_TAU_YEARS * 365.25

# ═══════════════════════════════════════════════════════════════════════════════
# 5. EVENT-COUNT GATING FRAMEWORK (Protocol Section 6.3)
# ═══════════════════════════════════════════════════════════════════════════════

#' Determine analytic tier based on observed death event count
#' @param n_events Integer: total death events (Outcome_OS == 1)
#' @return Named list with tier label, permitted analyses description,
#'         max_df, and boolean flags for what is permitted
get_analytic_tier <- function(n_events) {
  if (n_events < 80) {
    tier <- list(
      label        = "Inadequate (< 80 events)",
      tier_code    = 1,
      max_df       = 2,
      permit_mva   = FALSE,
      permit_rcs   = FALSE,
      permit_dca   = FALSE,
      permit_interaction = FALSE,
      permit_subgroup_inference = FALSE,
      description  = paste0(
        "KM curves + log-rank test only. Univariate Cox for LNR and N-stage. ",
        "No multivariable adjustment. Findings = hypothesis-generating."
      )
    )
  } else if (n_events < 130) {
    tier <- list(
      label        = "Minimum Viable (80-129 events)",
      tier_code    = 2,
      max_df       = 8,
      permit_mva   = TRUE,
      permit_rcs   = FALSE,   # Use linear terms, not RCS
      permit_dca   = TRUE,
      permit_interaction = FALSE,
      permit_subgroup_inference = FALSE,
      description  = paste0(
        "Restricted MVA: max 8 df, linear terms for LNR & TLE. ",
        "C-statistic + DCA + calibration. No interaction tests."
      )
    )
  } else if (n_events < 200) {
    tier <- list(
      label        = "Target (130-199 events)",
      tier_code    = 3,
      max_df       = 16,
      permit_mva   = TRUE,
      permit_rcs   = TRUE,
      permit_dca   = TRUE,
      permit_interaction = FALSE,
      permit_subgroup_inference = FALSE,
      description  = paste0(
        "Full MVA with RCS splines. Uno C-stat, DCA, calibration. ",
        "Fine-Gray for CSS. Subgroups descriptive only."
      )
    )
  } else {
    tier <- list(
      label        = "Full Power (>= 200 events)",
      tier_code    = 4,
      max_df       = 20,
      permit_mva   = TRUE,
      permit_rcs   = TRUE,
      permit_dca   = TRUE,
      permit_interaction = TRUE,
      permit_subgroup_inference = TRUE,
      description  = paste0(
        "All analyses + formal LNR x TLE interaction. ",
        "Subgroup analyses with inferential C-stat."
      )
    )
  }
  tier$n_events <- n_events
  return(tier)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 6. PUBLICATION THEME AND COLOURS
# ═══════════════════════════════════════════════════════════════════════════════

# -- Colour palettes --
COL_NSTAGE <- c("N1" = "#2166AC", "N2" = "#B2182B")

COL_LNR_TIERS <- c(
  "Low (<0.05)"              = "#4DAF4A",
  "Intermediate (0.05-0.20)" = "#FF7F00",
  "High (>0.20)"             = "#E41A1C"
)

COL_MODELS <- c(
  "Model A (N-stage)"          = "#2166AC",
  "Model B (LNR)"              = "#B2182B",
  "Model C2 (N-stage + TLE)"   = "#4393C3",
  "Model D2 (LNR + TLE)"       = "#D6604D",
  "Model C (N-stage Full MVA)"  = "#053061",
  "Model D (LNR Full MVA)"      = "#67001F"
)

COL_TLE <- c("Inadequate (<12)" = "#D95F02", "Adequate (>=12)" = "#1B9E77")

# -- Publication ggplot2 theme --
theme_publication <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      # Text
      plot.title       = element_text(face = "bold", size = base_size + 2, 
                                       hjust = 0),
      plot.subtitle    = element_text(size = base_size, colour = "grey40"),
      plot.caption     = element_text(size = base_size - 2, colour = "grey50",
                                       hjust = 0),
      # Axes
      axis.title       = element_text(face = "bold", size = base_size),
      axis.text        = element_text(size = base_size - 1),
      axis.line        = element_line(colour = "black", linewidth = 0.4),
      # Legend
      legend.title     = element_text(face = "bold", size = base_size - 1),
      legend.text      = element_text(size = base_size - 2),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key       = element_rect(fill = "white", colour = NA),
      legend.position  = "bottom",
      # Panel
      panel.border     = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      # Strip (facets)
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      # Margins
      plot.margin      = margin(10, 10, 10, 10)
    )
}

# Set as default theme
theme_set(theme_publication())

# -- Figure export settings --
FIG_WIDTH  <- 8    # inches
FIG_HEIGHT <- 6    # inches
FIG_DPI    <- 300  # dots per inch

# ═══════════════════════════════════════════════════════════════════════════════
# 7. HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

#' Save a ggplot figure to both PDF and PNG
#' @param plot ggplot object
#' @param filename Base filename without extension
#' @param width Width in inches (default FIG_WIDTH)
#' @param height Height in inches (default FIG_HEIGHT)
save_figure <- function(plot, filename, width = FIG_WIDTH, height = FIG_HEIGHT,
                        subdir = DIR_FIGURES) {
  ggsave(file.path(subdir, paste0(filename, ".pdf")), 
         plot = plot, width = width, height = height, dpi = FIG_DPI)
  ggsave(file.path(subdir, paste0(filename, ".png")), 
         plot = plot, width = width, height = height, dpi = FIG_DPI)
  cat("Saved:", filename, "(PDF + PNG)\n")
}

#' Save a table as CSV and return formatted kable
#' @param df Data frame
#' @param filename Base filename without extension
save_table <- function(df, filename, subdir = DIR_TABLES) {
  write.csv(df, file.path(subdir, paste0(filename, ".csv")), 
            row.names = FALSE)
  cat("Saved:", filename, ".csv\n")
}

#' Print a section header to console
section_header <- function(title) {
  width <- 78
  cat("\n", strrep("=", width), "\n", sep = "")
  cat("  ", title, "\n", sep = "")
  cat(strrep("=", width), "\n\n", sep = "")
}

#' Format p-value for display
format_pval <- function(p, digits = 3) {
  ifelse(p < 0.001, "< 0.001",
         ifelse(p < 0.01, sprintf("%.3f", p),
                sprintf(paste0("%.", digits, "f"), p)))
}

#' Wilson confidence interval for a proportion
wilson_ci <- function(x, n, alpha = ALPHA) {
  z <- qnorm(1 - alpha / 2)
  p_hat <- x / n
  denom <- 1 + z^2 / n
  centre <- (p_hat + z^2 / (2 * n)) / denom
  margin <- (z / denom) * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))
  c(lower = max(0, centre - margin), 
    estimate = p_hat, 
    upper = min(1, centre + margin))
}

cat("\n", strrep("=", 78), "\n", sep = "")
cat("  CONFIG LOADED SUCCESSFULLY\n")
cat("  Node-positive filter: ", NODE_POSITIVE_FILTER, "\n")
cat("  LNR cutpoints (Rosenberg): ", LNR_CUT_LOW_INTER, "/", LNR_CUT_INTER_HIGH, "\n")
cat("  TLE adequacy threshold: ", TLE_ADEQUATE_THRESHOLD, "\n")
cat("  Bootstrap replicates: ", N_BOOTSTRAP, "\n")
cat("  Random seed: ", SEED, "\n")
cat(strrep("=", 78), "\n")
