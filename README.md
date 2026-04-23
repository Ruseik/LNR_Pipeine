# Prognostic Comparison of LNR vs. N-Stage in Colorectal Cancer

## Overview
This repository contains a fully reproducible, modular R-based analytical pipeline designed to rigorously compare the prognostic performance of **Lymph Node Ratio (LNR)** and conventional **AJCC N-stage** in node-positive colorectal cancer (CRC). 

Conventional N-stage (N1/N2) is often limited by its dependence on the total number of lymph nodes examined (TLE). LNR—the ratio of metastatic to examined lymph nodes—offers a more nuanced metric that accounts for variability in lymph node yield. This pipeline provides a "plug-and-play" framework to evaluate these metrics across multiple statistical dimensions, including discrimination, calibration, and clinical utility.

---

## 🧬 Scientific Background
The prognostic adequacy of pathological nodal staging is a cornerstone of CRC management. However, inadequate lymph node yield (TLE < 12) can lead to "understaging." LNR has emerged as a potentially superior metric. 

**Key Analytical Pillars:**
- **Flexibility**: Modeling LNR as a continuous variable using Restricted Cubic Splines (RCS).
- **Comparison**: Three-layer modeling (Unadjusted, TLE-adjusted, and Fully adjusted).
- **Utility**: Assessment via Decision Curve Analysis (DCA) and Uno’s C-statistic.
- **Robustness**: Internal validation via bootstrapping and competing risks analysis (Fine-Gray).

---

## 🛠 Pipeline Architecture
The pipeline is implemented as a sequence of intuitively numbered R scripts. Each module is self-contained but inherits global settings from the master configuration.

| Script | Module | Description |
| :--- | :--- | :--- |
| `00_config.R` | **Configuration** | Global paths, packages, seed, and analytical constants (e.g., Rosenberg thresholds). |
| `01_data_prep.R` | **Data Preparation** | Cleaning, integrity checks, derivation of LNR tiers, and event-count gating. |
| `02_descriptive.R` | **Descriptive Stats** | Table 1 (SMD-based), TLE distribution, and LNR characterization. |
| `03_unadjusted.R` | **Unadjusted Survival** | Kaplan-Meier curves, Log-rank tests, and 1/3/5-year survival estimates. |
| `04_cox_models.R` | **Cox Regression** | Multi-layer model fitting (N-stage vs. LNR) with RCS integration. |
| `05_performance.R` | **Performance Metrics** | Uno's C-statistic, Calibration plots, and Decision Curve Analysis (DCA). |
| `06_reclass.R` | **Reclassification** | Cross-tabulation of N-stage vs. LNR tiers and RMST analysis. |
| `07_subgroups.R` | **Sensitivity Analysis** | Stratified analyses (TLE <12 vs ≥12) and alternative LNR cut-points. |
| `08_css_analysis.R`| **Competing Risks** | Fine-Gray models for Cancer-Specific Survival (CSS). |
| `09_reporting.R` | **Final Outputs** | Compilation of publication-ready tables and multi-panel figures. |

---

## 🚀 Key Features

### 1. Event-Count Gating Framework
To ensure statistical rigor, the pipeline automatically scales model complexity based on the observed death events:
- **Inadequate (<80 events)**: KM and Univariate only.
- **Minimum Viable (80-129 events)**: Restricted MVA (max 8 df).
- **Target (130-199 events)**: Full MVA with RCS splines.
- **Full Power (≥200 events)**: Includes formal LNR × TLE interaction tests.

### 2. Advanced Survival Modeling
- **Non-linearity**: Uses `rms::rcs` for Restricted Cubic Splines to capture non-linear effects of LNR and TLE.
- **Clinical Utility**: `dcurves::dca` for Decision Curve Analysis to determine net benefit across threshold probabilities.
- **Validation**: Bootstrap-based optimism-corrected C-statistics.

---

## 📋 Data Requirements
The pipeline expects a CSV dataset (default: `CRCALVS_8_1.csv`) with the following key variables:
- `Index`: Patient identifier.
- `Time_survival`: Follow-up time (days).
- `Outcome_OS`: Overall survival status (1 = death, 0 = censored).
- `LN_total`: Total lymph nodes examined (TLE).
- `LN_positive`: Number of metastatic lymph nodes.
- `AG`, `GN`: Age and Gender.
- `ST_main`, `SN_main`: Pathological T and N stages.
- `Grade`, `LVI`: Tumor grade and Lymphovascular invasion.

---

## ⚙️ Setup & Usage

### Prerequisites
- **R version**: 4.0.0+ recommended.
- **Core Packages**: `survival`, `rms`, `survminer`, `cmprsk`, `dcurves`, `mice`, `tableone`, `patchwork`.

### Running the Pipeline
1.  **Configure Paths**: Open `Analysis/00_config.R` and update `DIR_PROJECT` and `FILE_RAW_DATA` to match your local environment.
2.  **Execute in Sequence**:
    ```R
    # Example: Running the full pipeline from R console
    source("Analysis/00_config.R")
    source("Analysis/01_data_preparation.R")
    # ... continue through 09_tables_figures.R
    ```

---

## 📂 Project Structure
```text
Analysis/
├── 00_config.R              # Master settings
├── 01_data_preparation.R    # Data cleaning & gating
├── ...                      # Modular analysis scripts
├── Outputs/                 # Auto-generated results
│   ├── Figures/             # KM plots, DCA curves, Calibration plots
│   ├── Tables/              # Table 1, Model comparison, C-stats
│   ├── Models/              # Saved RData model objects
│   └── Data/                # Cleaned analytic datasets
└── Methodology/             # Pre-specified methodology docs
```

---

## 📝 Abstract
The prognostic adequacy of conventional pathological nodal staging (N-stage) in node-positive colorectal cancer is limited by its inability to account for variability in lymph node yield. The lymph node ratio (LNR), defined as the proportion of metastatic to examined lymph nodes, has emerged as a potentially superior metric. We developed a fully reproducible, modular R-based analytical pipeline implemented on intuitively numbered scripts to rigorously compare the prognostic performance of N-stage and LNR across multiple statistical dimensions. This plug-and-play framework requires only formatted dataset substitution, with all analyses automated. 

Findings support LNR as a clinically informative alternative to N-stage within a scalable, transferable analytic framework.

---
*Developed for robust comparative inference in surgical oncology.*
