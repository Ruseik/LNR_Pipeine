# Implementation Plan: LNR vs N-Stage Prognostic Analysis — R Scripts

## Goal

Translate the Version 2.0 Research Protocol into a modular set of R scripts that execute the full analytic pipeline: data preparation → descriptive statistics → unadjusted survival → three-layer model comparison → performance metrics → reclassification → subgroup/sensitivity analyses. Primary outcome: **Overall Survival (OS)**. Secondary outcome: **Cancer-Specific Survival (CSS)**. DFS is excluded per user instruction.

---

## User Review Required

> [!IMPORTANT]
> **Dataset–Protocol Variable Mismatch**
> The protocol specifies several covariates that are **NOT present** in the dataset [CRCALVS_8_1.csv](file:///d:/Google_SSD_RAM/_OneDrive/Research/Project_ColoRectal/Analysis_Project_ColoRectal/CRCALVS_8_1.csv):
>
> | Protocol Covariate | Status in Dataset |
> |---|---|
> | Perineural invasion (PNI) | ❌ Not available |
> | Resection margin (R0/R1) | ❌ Not available |
> | Tumour site (colon vs rectum) | ❌ Not available |
> | Calendar period / year of surgery | ❌ Not available |
> | Chemotherapy receipt/timing | ❌ Not available |
> | Centre/site | ❌ Not available (assumed single-centre) |
>
> **Consequence**: The full MVA model (Layer 3) will use the **available** covariates only: Age, Sex, T-stage, Grade, LVI, and TLE. The model degrees of freedom are reduced, which lowers the event-count requirement. This must be acknowledged as a limitation.

> [!WARNING]
> **Small Sample Size**
> The dataset has **138 total records**. After filtering to node-positive only (SN_main ≥ 1), the analytic cohort will be substantially smaller (~40-55 patients based on inspection). The death event count is likely to fall in the **Inadequate** or **Minimum Viable** tier of the event-count gating framework. The plan includes the gating logic to automatically scale model complexity.

> [!IMPORTANT]
> **LNR Encoding**
> LNR in the dataset is stored as a **percentage (0–100)**, not a proportion (0–1). The scripts will divide by 100 to convert to the 0–1 scale used in the protocol and the Rosenberg cutpoints (0.05/0.20).

> [!IMPORTANT]
> **N-Stage Coding**
> `SN_main` is coded as 0/1/2, mapping to N0/N1/N2. The AJCC sub-categories (N1a/N1b/N1c/N2a/N2b) are **not available** — analysis uses the three-tier N0/N1/N2 classification (N0 filtered out). No N1c identification is possible, so the N1c sensitivity analysis is not feasible.

---

## Proposed Changes

### Script Architecture

All scripts are written to `d:\Google_SSD_RAM\_OneDrive\Research\Project_ColoRectal\LNR_20260328\Analysis\`. A master config script is sourced by all subsequent scripts.

---

### Script 00: Master Configuration

#### [NEW] [00_config.R](file:///d:/Google_SSD_RAM/_OneDrive/Research/Project_ColoRectal/LNR_20260328/Analysis/00_config.R)

- Package loading: `survival`, `rms`, `cmprsk`, `dcurves`, `mice`, `ggplot2`, `survminer`, `pROC`, `boot`, `tableone`, `knitr`, `patchwork`
- Global paths: data file, output directory
- Colour palettes & theme settings for publication-quality figures
- LNR tier cutpoints: `c(0.05, 0.20)` (Rosenberg)
- TLE adequacy threshold: 12
- Event-count gating tier definitions and a function to return the permitted tier
- Random seed for reproducibility: `set.seed(2024)`
- Uno C-stat truncation quantile: 0.75
- Bootstrap replicates: 1000
- DCA threshold range: `seq(0.10, 0.60, by = 0.01)`

---

### Script 01: Data Preparation

#### [NEW] [01_data_preparation.R](file:///d:/Google_SSD_RAM/_OneDrive/Research/Project_ColoRectal/LNR_20260328/Analysis/01_data_preparation.R)

**Steps:**
1. Read CSV, drop useless column `X`
2. **Convert LNR from percentage to proportion** (`LNR / 100`)
3. **Validate data integrity**: check `LN_positive ≤ LN_total`, flag any violations
4. **Recompute LNR** from `LN_positive / LN_total` and compare with stored LNR — report discrepancies
5. **Filter to node-positive cohort**: `SN_main >= 1` (i.e., N1 or N2)
6. **Create derived variables**:
   - `LNR_prop` — recomputed LNR as proportion (0–1)
   - `LNR_tier` — categorical: Low (< 0.05), Intermediate (0.05–0.20), High (> 0.20) per Rosenberg
   - `N_stage` — factor from `SN_main`: "N1" (1), "N2" (2)
   - `T_stage_group` — per protocol: pT1/pT2 vs pT3 vs pT4
   - `TLE_adequate` — binary: TLE ≥ 12 vs < 12
   - `Grade_f` — factor: 0 = "Low", 1 = "High"
   - `LVI_f` — factor: 0 = "Absent", 1 = "Present"
   - `Sex_f` — factor from `GN`: 0 = "Male", 1 = "Female"
7. **Convert survival time from days to months** (÷ 30.44) for display, keep days for modelling
8. **Assess missingness**: report proportion missing per variable
9. **Apply event-count gating**: count deaths (`Outcome_OS == 1`) in node-positive cohort → assign analytic tier
10. **Save** cleaned analytic dataset as `.RData` and [.csv](file:///d:/Google_SSD_RAM/_OneDrive/Research/Project_ColoRectal/Analysis_Project_ColoRectal/CRCALVS_8_1.csv)
11. **Generate CONSORT-style flow text**: total records → exclusions → final analytic N

---

### Script 02: Phase 1 — Descriptive Statistics (SO6)

#### [NEW] [02_descriptive.R](file:///d:/Google_SSD_RAM/_OneDrive/Research/Project_ColoRectal/LNR_20260328/Analysis/02_descriptive.R)

**Steps:**
1. **Table 1**: Baseline characteristics by LNR tier and by N-stage, using `tableone::CreateTableOne` with SMD instead of p-values
2. **TLE characterisation**:
   - Distribution: median, IQR, range, histogram
   - Proportion with TLE < 12 and TLE < 8 with 95% Wilson CI
   - (Calendar period trend not possible — no date variable)
3. **LNR distribution**:
   - Histogram of continuous LNR
   - Distribution across Rosenberg tiers: counts, percentages, CI
   - Comparison of LNR tier quantiles against published Western distributions
4. **Outcome summary**: number of deaths (OS), DSS events, median follow-up (reverse KM)
5. **Patient flow diagram** (text-based CONSORT)

---

### Script 03: Phase 2 — Unadjusted Survival (SO1)

#### [NEW] [03_unadjusted_survival.R](file:///d:/Google_SSD_RAM/_OneDrive/Research/Project_ColoRectal/LNR_20260328/Analysis/03_unadjusted_survival.R)

**Steps:**
1. **KM curves — N-stage** (N1 vs N2): plot + log-rank test + 1/3/5-year survival estimates with 95% CI
2. **KM curves — LNR tiers** (Low/Intermediate/High): plot + log-rank + 1/3/5-year estimates
3. **RMST at 5 years** for each group (N-stage, LNR tier)
4. **Continuous LNR dose-response**: martingale residual plot from Cox null model to assess linearity
5. **Univariate Cox models**: LNR continuous, LNR tiers, N-stage — HRs with 95% CI

---

### Script 04: Phase 3 — Three-Layer Cox Model Comparison (SO2, SO3)

#### [NEW] [04_cox_models.R](file:///d:/Google_SSD_RAM/_OneDrive/Research/Project_ColoRectal/LNR_20260328/Analysis/04_cox_models.R)

**Steps:**

Implements the six-model three-layer symmetric comparison:

| Layer | Model | Formula |
|---|---|---|
| 1 | A | `Surv(Time_survival, Outcome_OS) ~ N_stage` |
| 1 | B | `Surv(Time_survival, Outcome_OS) ~ rcs(LNR_prop, 3)` |
| 2 | C₂ | `~ N_stage + rcs(LN_total, 3)` |
| 2 | D₂ | `~ rcs(LNR_prop, 3) + rcs(LN_total, 3)` |
| 3 | C | `~ N_stage + rcs(LN_total, 3) + AG + Sex_f + T_stage_group + Grade_f + LVI_f` |
| 3 | D | `~ rcs(LNR_prop, 3) + rcs(LN_total, 3) + AG + Sex_f + T_stage_group + Grade_f + LVI_f` |

**Event-count gating logic** (automatic):
- If deaths < 80: only Models A and B (univariate); report as hypothesis-generating
- If deaths 80–129: reduce to max 8 df — use linear terms for LNR and TLE, not RCS
- If deaths ≥ 130: full RCS specification as above
- If deaths ≥ 200: add LNR × TLE_adequate interaction model

**For each model:**
- Cox PH fit via `rms::cph()` for RCS and validation compatibility
- Report coefficients, HR, 95% CI, p-values
- PH assumption testing via `cox.zph()` — if violated, apply stratification or time-interaction

**Association plot**: log-HR as function of continuous LNR with 95% CI bands (from RCS model)

---

### Script 05: Phase 3 — Performance Metrics (SO2)

#### [NEW] [05_performance_metrics.R](file:///d:/Google_SSD_RAM/_OneDrive/Research/Project_ColoRectal/LNR_20260328/Analysis/05_performance_metrics.R)

**Primary co-endpoints:**
1. **Uno's C-statistic** for all 6 models, truncated at 75th percentile of event times
   - ΔC between Model D and Model C (Layer 3) with bootstrap 95% CI (1000 resamples, BCa)
   - Pre-specified decision threshold: ΔC ≥ 0.05
2. **Decision-Curve Analysis (DCA)** using `dcurves::dca()`:
   - Models C vs D (Layer 3) vs treat-all vs treat-none
   - Threshold probability range: 10–60%
   - Net benefit comparison plot

**Primary secondary endpoints:**
3. **Calibration** at 3 and 5 years:
   - Calibration-in-the-large
   - Calibration slope
   - Calibration plot (predicted vs observed by deciles)
4. **AIC and BIC** for all 6 models

**Exploratory (appendix):**
5. Continuous NRI (Pencina 2011) — appendix only
6. IDI — appendix only

**Internal validation:**
7. **Bootstrap optimism-corrected C-statistic** using `rms::validate()` (.632+ method, 1000 resamples)

---

### Script 06: Phase 4 — Reclassification Analysis (SO4)

#### [NEW] [06_reclassification.R](file:///d:/Google_SSD_RAM/_OneDrive/Research/Project_ColoRectal/LNR_20260328/Analysis/06_reclassification.R)

1. **Cross-tabulation**: N-stage (N1/N2) × LNR tier (Low/Intermediate/High) — counts and proportions
2. **Reclassification proportion**: patients whose LNR tier implies a different risk than their N-stage
3. **Survival by reclassification cell**: KM curves for each cell of the cross-tabulation
4. **RMST comparison**: do reclassified patients cluster with LNR-assigned or N-stage group?

---

### Script 07: Phase 5 — Subgroup & Sensitivity Analyses (SO5)

#### [NEW] [07_subgroup_sensitivity.R](file:///d:/Google_SSD_RAM/_OneDrive/Research/Project_ColoRectal/LNR_20260328/Analysis/07_subgroup_sensitivity.R)

**Pre-specified subgroups** (descriptive unless ≥ 200 events):
1. TLE < 12 vs ≥ 12: C-statistic for LNR vs N-stage within each stratum
2. T3/T4 vs T1/T2
3. N1 vs N2

**Pre-specified sensitivity analyses:**
1. Complete case analysis (no imputation needed if missingness < 5%)
2. Landmark analysis at 6 months: exclude events in first 180 days
3. Binary TLE (< 12 vs ≥ 12) replacing RCS TLE in MVA
4. Alternative LNR thresholds: 0.05/0.40 and 0.20 binary cutpoint

---

### Script 08: CSS Analysis (Secondary Outcome)

#### [NEW] [08_css_analysis.R](file:///d:/Google_SSD_RAM/_OneDrive/Research/Project_ColoRectal/LNR_20260328/Analysis/08_css_analysis.R)

1. **Check CSS data completeness**: proportion of deceased with `Outcome_DSS` available
   - If ≥ 85% cause-of-death data: proceed with Fine-Gray competing risks model
   - If < 85%: descriptive analysis only
2. **Fine-Gray subdistribution hazard models** for CSS:
   - LNR tier vs N-stage as competing-risks predictors
   - Non-cancer death as competing event
3. **Cumulative incidence curves** by LNR tier and N-stage

---

### Script 09: Tables & Figures Compilation

#### [NEW] [09_tables_figures.R](file:///d:/Google_SSD_RAM/_OneDrive/Research/Project_ColoRectal/LNR_20260328/Analysis/09_tables_figures.R)

1. Compile all publication-ready tables (Table 1, model comparison table, C-statistic table, reclassification table)
2. Multi-panel figure assembly using `patchwork`
3. Export tables as CSV/HTML and figures as PDF/PNG (300 dpi)

---

## Key Design Decisions

### Degrees of Freedom Budget (Adjusted for Available Variables)

With the reduced covariate set (no PNI, margin, site, calendar period):

| Variable | df |
|---|---|
| LNR (RCS 3 knots) or N-stage (1 level) | 2 or 1 |
| TLE (RCS 3 knots) | 2 |
| Age (linear) | 1 |
| Sex | 1 |
| T-stage group (3 levels) | 2 |
| Grade | 1 |
| LVI | 1 |
| **Total (LNR model)** | **10** |
| **Total (N-stage model)** | **9** |

At EPV = 10, this requires **~100 death events** — likely achievable depending on the node-positive subset size and mortality rate.

### Multiple Imputation

- Applied via `mice` with m ≥ 20 only if any covariate has > 5% missingness
- If all covariates are complete (likely for this small clean dataset), complete case is the primary analysis

### Automatic Gating

All scripts query the event count from `01_data_preparation.R` output and automatically scale complexity. If events < 80, scripts 04/05/06/07 reduce to univariable analyses and clearly label outputs as hypothesis-generating.

---

## Verification Plan

### Automated Tests

Since this is an R statistical analysis (not software development), verification is through execution and output inspection:

1. **Run all scripts in sequence**:
   ```
   Rscript 00_config.R
   Rscript 01_data_preparation.R
   Rscript 02_descriptive.R
   Rscript 03_unadjusted_survival.R
   Rscript 04_cox_models.R
   Rscript 05_performance_metrics.R
   Rscript 06_reclassification.R
   Rscript 07_subgroup_sensitivity.R
   Rscript 08_css_analysis.R
   Rscript 09_tables_figures.R
   ```
2. **Check**: no errors/warnings, all output files generated in `Outputs/` directory
3. **Validate**: LNR recomputation matches stored values, node-positive filter correct, CONSORT flow numbers add up

### Manual Verification

- Review KM curves for clinical plausibility (N2 should have worse survival than N1; High LNR worse than Low)
- Review C-statistic values are in plausible range (0.5–0.9)
- Review calibration plots for obvious miscalibration
- Cross-check reclassification table: row/column totals match analytic N
