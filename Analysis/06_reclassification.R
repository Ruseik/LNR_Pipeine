################################################################################
#  06_reclassification.R — Phase 4: Reclassification Analysis
#  Protocol Section 7.5 | SO4
#
#  Cross-tabulation of N-stage × LNR tier, reclassification proportion,
#  and survival validation of reclassified patients
################################################################################

source(file.path("d:/Google_SSD_RAM/_OneDrive/Research",
                 "Project_ColoRectal/LNR_20260328/Analysis/00_config.R"))

load(file.path(DIR_DATA, "analytic_dataset.RData"))

section_header("SCRIPT 06: RECLASSIFICATION ANALYSIS (Phase 4)")

cat("Analytic cohort: N =", nrow(df), "| OS events:", 
    sum(df$Outcome_OS == 1), "\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 1. CROSS-TABULATION: N-STAGE × LNR TIER
# ═══════════════════════════════════════════════════════════════════════════════

section_header("N-STAGE × LNR TIER CROSS-TABULATION")

cross_tab <- table(df$N_stage, df$LNR_tier)
cat("Counts:\n")
print(addmargins(cross_tab))

cat("\nRow percentages (within each N-stage):\n")
row_pct <- round(prop.table(cross_tab, margin = 1) * 100, 1)
print(row_pct)

cat("\nColumn percentages (within each LNR tier):\n")
col_pct <- round(prop.table(cross_tab, margin = 2) * 100, 1)
print(col_pct)

# ═══════════════════════════════════════════════════════════════════════════════
# 2. RECLASSIFICATION ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

section_header("RECLASSIFICATION ANALYSIS")

# Define "concordant" vs "discordant" classification
# N1 patients → LOW/INTERMEDIATE LNR = concordant (expected lower risk)
#             → HIGH LNR = reclassified upward (higher risk than N-stage implies)
# N2 patients → HIGH LNR = concordant (expected higher risk)
#             → LOW/INTERMEDIATE LNR = reclassified downward (lower risk)

df$reclass_group <- NA_character_

# N1 patients
n1_mask <- df$N_stage == "N1"
df$reclass_group[n1_mask & df$LNR_tier == "Low (<0.05)"] <- "N1_LNR-Low (concordant)"
df$reclass_group[n1_mask & df$LNR_tier == "Intermediate (0.05-0.20)"] <- "N1_LNR-Intermediate"
df$reclass_group[n1_mask & df$LNR_tier == "High (>0.20)"] <- "N1_LNR-High (upclassified)"

# N2 patients
n2_mask <- df$N_stage == "N2"
df$reclass_group[n2_mask & df$LNR_tier == "Low (<0.05)"] <- "N2_LNR-Low (downclassified)"
df$reclass_group[n2_mask & df$LNR_tier == "Intermediate (0.05-0.20)"] <- "N2_LNR-Intermediate (downclassified)"
df$reclass_group[n2_mask & df$LNR_tier == "High (>0.20)"] <- "N2_LNR-High (concordant)"

cat("Reclassification groups:\n")
print(table(df$reclass_group, useNA = "ifany"))

# Proportions
n_total <- nrow(df)
n_concordant <- sum(grepl("concordant", df$reclass_group, ignore.case = TRUE), 
                    na.rm = TRUE)
n_reclassified <- n_total - n_concordant

cat("\nConcordant (N-stage and LNR agree on risk level):", n_concordant, 
    "(", round(100 * n_concordant / n_total, 1), "%)\n")
cat("Reclassified (N-stage and LNR disagree):", n_reclassified,
    "(", round(100 * n_reclassified / n_total, 1), "%)\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 3. SURVIVAL BY RECLASSIFICATION CELL
# ═══════════════════════════════════════════════════════════════════════════════

section_header("SURVIVAL BY RECLASSIFICATION CELL")

# Create combined cell variable
df$reclass_cell <- paste0(df$N_stage, " / ", df$LNR_tier)
df$reclass_cell <- factor(df$reclass_cell)

cat("Cell counts and events:\n")
cell_summary <- df %>%
  group_by(reclass_cell) %>%
  summarise(
    N = n(),
    Events_OS = sum(Outcome_OS == 1),
    Median_surv_days = median(Time_survival),
    .groups = "drop"
  ) %>%
  as.data.frame()
print(cell_summary)
save_table(cell_summary, "reclassification_cell_summary")

# KM curves by reclassification cell (only for cells with sufficient data)
cells_with_data <- cell_summary$reclass_cell[cell_summary$N >= 3 & 
                                               cell_summary$Events_OS >= 1]

if (length(cells_with_data) >= 2) {
  df_km_reclass <- df[df$reclass_cell %in% cells_with_data, ]
  
  km_reclass <- survfit(Surv(Time_survival, Outcome_OS) ~ reclass_cell, 
                         data = df_km_reclass)
  
  # Custom colours
  n_cells <- length(cells_with_data)
  cell_colours <- viridis(n_cells)
  
  p_km_reclass <- ggsurvplot(
    km_reclass, data = df_km_reclass,
    palette = cell_colours,
    risk.table = TRUE,
    risk.table.height = 0.30,
    pval = TRUE,
    conf.int = FALSE,
    xlab = "Time (days)",
    ylab = "Overall Survival Probability",
    title = "Overall Survival by N-Stage × LNR Tier Classification",
    subtitle = "Each cell represents a unique concordance/discordance group",
    legend.title = "N-Stage / LNR Tier",
    ggtheme = theme_publication(),
    break.time.by = 365.25,
    legend = "right"
  )
  
  pdf(file.path(DIR_FIGURES, "Fig_KM_reclassification_cells.pdf"),
      width = FIG_WIDTH + 2, height = FIG_HEIGHT + 2)
  print(p_km_reclass)
  dev.off()
  png(file.path(DIR_FIGURES, "Fig_KM_reclassification_cells.png"),
      width = FIG_WIDTH + 2, height = FIG_HEIGHT + 2, units = "in", res = FIG_DPI)
  print(p_km_reclass)
  dev.off()
  cat("Saved: Fig_KM_reclassification_cells\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 4. RMST BY RECLASSIFICATION CELL
# ═══════════════════════════════════════════════════════════════════════════════

section_header("RMST BY RECLASSIFICATION CELL")

tau <- min(RMST_TAU_DAYS, max(df$Time_survival, na.rm = TRUE))

rmst_cells <- data.frame()
for (cell in levels(df$reclass_cell)) {
  sub <- df[df$reclass_cell == cell, ]
  if (nrow(sub) >= 2 && sum(sub$Outcome_OS == 1) >= 1) {
    km_sub <- survfit(Surv(Time_survival, Outcome_OS) ~ 1, data = sub)
    summ <- summary(km_sub, rmean = tau)$table
    # Use grep to find rmean fields (names vary across survival pkg versions)
    rmean_idx <- grep("rmean", names(summ))[1]
    se_idx    <- grep("se.*rmean", names(summ))[1]
    rmst_cells <- rbind(rmst_cells, data.frame(
      Cell = cell,
      N = nrow(sub),
      Events = sum(sub$Outcome_OS == 1),
      RMST_days = round(unname(summ[rmean_idx]), 0),
      RMST_years = round(unname(summ[rmean_idx]) / 365.25, 2),
      SE = round(unname(summ[se_idx]), 0),
      stringsAsFactors = FALSE
    ))
  }
}

cat("RMST at", round(tau / 365.25, 1), "years by reclassification cell:\n")
print(rmst_cells)
save_table(rmst_cells, "RMST_reclassification_cells")

# ═══════════════════════════════════════════════════════════════════════════════
# 5. MOSAIC/HEATMAP: RECLASSIFICATION WITH SURVIVAL
# ═══════════════════════════════════════════════════════════════════════════════

section_header("RECLASSIFICATION HEATMAP")

# Create a heatmap showing cell counts with survival overlay
cross_df <- as.data.frame(cross_tab)
names(cross_df) <- c("N_Stage", "LNR_Tier", "Count")

# Calculate survival rate per cell
cross_df$Mortality_rate <- NA
for (i in seq_len(nrow(cross_df))) {
  sub <- df[df$N_stage == cross_df$N_Stage[i] & 
             df$LNR_tier == cross_df$LNR_Tier[i], ]
  if (nrow(sub) > 0) {
    cross_df$Mortality_rate[i] <- round(100 * sum(sub$Outcome_OS == 1) / nrow(sub), 1)
  }
}

p_heatmap <- ggplot(cross_df, aes(x = LNR_Tier, y = N_Stage, fill = Mortality_rate)) +
  geom_tile(colour = "white", linewidth = 1.5) +
  geom_text(aes(label = paste0("N = ", Count, "\n", 
                                ifelse(is.na(Mortality_rate), "", 
                                       paste0(Mortality_rate, "% died")))),
            colour = "white", fontface = "bold", size = 4) +
  scale_fill_gradient(low = "#2166AC", high = "#B2182B", na.value = "grey80",
                      name = "Mortality\nRate (%)") +
  labs(title = "Reclassification Heatmap: N-Stage × LNR Tier",
       subtitle = "Cell size and mortality rate for OS",
       x = "LNR Tier (Rosenberg)", y = "N-Stage") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
save_figure(p_heatmap, "Fig_reclassification_heatmap")

# ═══════════════════════════════════════════════════════════════════════════════
# 6. SAVE
# ═══════════════════════════════════════════════════════════════════════════════

# Full reclassification table for manuscript
# Create matching composite key for merge
cross_df$Cell <- paste0(cross_df$N_Stage, " / ", cross_df$LNR_Tier)
reclass_table <- merge(cross_df, rmst_cells, by = "Cell", all.x = TRUE)
cross_df$Cell <- NULL   # remove temp key
save_table(cross_df, "reclassification_table_full")

save(cross_tab, cell_summary, rmst_cells,
     file = file.path(DIR_MODELS, "phase4_reclassification.RData"))

cat("\n", strrep("=", 78), "\n")
cat("  SCRIPT 06 COMPLETE\n")
cat(strrep("=", 78), "\n")
