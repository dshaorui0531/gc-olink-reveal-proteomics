suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

root_dir <- Sys.getenv("GC_OLINK_PROJECT_DIR", ".")
table_dir <- file.path(root_dir, "results", "tables")
out_dir <- file.path(root_dir, "results", "individual_panels_pub")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path(root_dir, ".fontconfig-cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = cache_dir)

figure_font <- "Helvetica"
theme_set(
  theme_classic(base_size = 6.8, base_family = figure_font) +
    theme(
      text = element_text(family = figure_font, colour = "black", lineheight = 0.95),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(colour = "black"),
      legend.title = element_text(size = 6.2),
      legend.text = element_text(size = 5.9, family = figure_font, lineheight = 0.95, margin = margin(r = 2)),
      legend.key.size = unit(3.0, "mm"),
      legend.spacing.x = unit(1.0, "mm"),
      legend.spacing.y = unit(1.0, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = 6.4, face = "bold"),
      plot.title = element_text(size = 7.2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 6.2, colour = "grey30"),
      panel.grid.major.y = element_line(linewidth = 0.18, colour = "grey90"),
      panel.grid.major.x = element_line(linewidth = 0.18, colour = "grey92")
    )
)

save_pub <- function(plot, stem, width_mm = 88, height_mm = 70, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svg_path <- file.path(out_dir, paste0(stem, ".svg"))
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  tiff_path <- file.path(out_dir, paste0(stem, ".tiff"))
  svglite(svg_path, width = w, height = h, bg = "white", fix_text_size = FALSE, system_fonts = list(sans = figure_font, serif = figure_font, mono = figure_font))
  print(plot)
  dev.off()
  pdf(pdf_path, width = w, height = h, family = figure_font, useDingbats = FALSE)
  print(plot)
  dev.off()
  agg_tiff(tiff_path, width = w, height = h, units = "in", res = dpi, background = "white")
  print(plot)
  dev.off()
}

marker_auc <- read.csv(file.path(table_dir, "traditional_marker_auc_E_vs_L.csv"), check.names = FALSE)
candidate_robust <- read.csv(file.path(table_dir, "candidate_clinical_adjustment_robustness.csv"), check.names = FALSE)
progression_robust <- read.csv(file.path(table_dir, "progression_candidate_multivariable_robustness.csv"), check.names = FALSE)
calibration <- read.csv(file.path(table_dir, "model_calibration_summary.csv"), check.names = FALSE)
dca <- read.csv(file.path(table_dir, "decision_curve_analysis_data.csv"), check.names = FALSE)

p_marker <- marker_auc |>
  mutate(
    Marker = factor(Marker, levels = Marker[order(AUC)]),
    label = paste0("n=", N)
  ) |>
  ggplot(aes(AUC, Marker)) +
  geom_vline(xintercept = 0.5, linewidth = 0.35, linetype = "dashed", colour = "grey45") +
  geom_errorbar(aes(xmin = AUC_CI_low, xmax = AUC_CI_high), orientation = "y", width = 0, linewidth = 0.45, colour = "grey40") +
  geom_point(size = 2.2, colour = "#2A6FBB") +
  geom_text(aes(label = label), nudge_x = 0.035, size = 2.1, hjust = 0, family = figure_font) +
  coord_cartesian(xlim = c(0.25, 1.05), clip = "off") +
  scale_x_continuous(breaks = seq(0.4, 1.0, 0.2)) +
labs(title = "Patient-only tumor marker AUCs", subtitle = "T3-T4 vs T1-T2 gastric cancer", x = "AUC (95% CI)", y = NULL) +
  theme(plot.margin = margin(5.5, 18, 5.5, 5.5))

save_pub(p_marker, "Fig4a_patient_only_marker_auc", 88, 62)

p_age_sex <- candidate_robust |>
  mutate(
    Comparison = recode(
      Comparison,
      E_vs_H = "T1-T2 GC vs healthy",
      L_vs_E = "T3-T4 vs T1-T2 GC",
      L_vs_H = "T3-T4 GC vs healthy"
    ),
    Retained = ifelse(Retained_age_sex_FDR_0.05, "FDR < 0.05 after adjustment", "Direction only")
  ) |>
  ggplot(aes(unadjusted_logFC, age_sex_logFC)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey45") +
  geom_point(aes(fill = Retained), shape = 21, size = 1.8, stroke = 0.25, colour = "black", alpha = 0.9) +
  facet_wrap(~ Comparison, nrow = 1) +
  scale_fill_manual(values = c("FDR < 0.05 after adjustment" = "#2A6FBB", "Direction only" = "#D8D8D8")) +
  labs(title = "Age- and sex-adjusted robustness", x = "Unadjusted log2 fold change", y = "Adjusted log2 fold change", fill = NULL) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 5.5, 9, 9)
  )

save_pub(p_age_sex, "Fig4b_age_sex_candidate_robustness", 178, 72)

signature_proteins <- c("MUC16", "NXPH3", "NPY", "ITGA11", "DSG3")

p_multivar <- progression_robust |>
  filter(Protein %in% signature_proteins) |>
  select(Protein, unadjusted_FDR, age_sex_FDR, demographic_FDR, full_clinical_FDR) |>
  pivot_longer(
    cols = ends_with("_FDR"),
    names_to = "Model",
    values_to = "FDR"
  ) |>
  mutate(
    Protein = factor(Protein, levels = rev(signature_proteins)),
    Model = recode(
      Model,
      unadjusted_FDR = "Unadjusted",
      age_sex_FDR = "Age + sex",
      demographic_FDR = "Demographic/lifestyle",
      full_clinical_FDR = "Full clinical"
    ),
    Model = factor(Model, levels = c("Unadjusted", "Age + sex", "Demographic/lifestyle", "Full clinical")),
    neglog10 = -log10(pmax(FDR, 1e-12)),
    Significant = FDR < 0.05
  ) |>
  ggplot(aes(neglog10, Protein, colour = Model)) +
  geom_vline(xintercept = -log10(0.05), linewidth = 0.35, linetype = "dashed", colour = "grey45") +
  geom_point(aes(shape = Significant), size = 2.1, position = position_dodge(width = 0.45)) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "FDR >= 0.05")) +
  scale_colour_manual(values = c(
    "Unadjusted" = "#111111",
    "Age + sex" = "#2A6FBB",
    "Demographic/lifestyle" = "#33A6A6",
    "Full clinical" = "#C45A4A"
  )) +
  labs(title = "Invasion-depth candidates under multivariable adjustment", x = expression(-log[10]("FDR")), y = NULL, colour = NULL, shape = NULL) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.width = unit(5.5, "mm"),
    plot.margin = margin(5.5, 5.5, 12, 5.5)
  ) +
  guides(
    colour = guide_legend(order = 1, nrow = 1, byrow = TRUE),
    shape = guide_legend(order = 2, nrow = 1, byrow = TRUE)
  )

save_pub(p_multivar, "Fig4c_progression_multivariable_robustness", 124, 86)

dca_keep <- dca |>
  filter(Threshold <= 0.85) |>
  mutate(
    Model = recode(
      Model,
      "Olink early signature" = "Olink T1-T2",
      "Clinical covariates: age + sex" = "Age+sex",
      "Olink progression signature" = "Olink inv.-depth",
      "Clinical covariates: age + sex + BMI" = "Age+sex+BMI"
    ),
    Model = factor(Model, levels = c("Olink T1-T2", "Age+sex", "Olink inv.-depth", "Age+sex+BMI")),
    Comparison = recode(Comparison, E_vs_H = "T1-T2 GC vs healthy", L_vs_E = "T3-T4 vs T1-T2 GC")
  )

p_dca <- ggplot(dca_keep, aes(Threshold, Net_benefit, colour = Model)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey65") +
  geom_line(linewidth = 0.55) +
  facet_wrap(~ Comparison, nrow = 1) +
  scale_colour_manual(values = c(
    "Olink T1-T2" = "#2A6FBB",
    "Age+sex" = "#9A9A9A",
    "Olink inv.-depth" = "#C45A4A",
    "Age+sex+BMI" = "#B8B8B8"
  ), na.value = "grey50") +
  coord_cartesian(xlim = c(0.05, 0.88), clip = "off") +
  labs(title = "Decision-curve summary", x = "Threshold probability", y = "Net benefit", colour = NULL) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(7, "mm"),
    plot.margin = margin(5.5, 6, 5.5, 5.5)
  ) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

p_brier <- calibration |>
  mutate(
    Model = recode(
      Model,
      "Olink early signature" = "Olink T1-T2",
      "Olink progression signature" = "Olink inv.-depth",
      "Clinical covariates: age + sex" = "Age+sex",
      "Clinical covariates: age + sex + BMI" = "Age+sex+BMI"
    ),
    Model = factor(Model, levels = Model[order(Brier_score)])
  ) |>
  ggplot(aes(Brier_score, Model)) +
  geom_col(width = 0.62, fill = "#6F8FB7") +
  geom_text(aes(label = sprintf("%.3f", Brier_score)), hjust = -0.08, size = 2.1, family = figure_font) +
  coord_cartesian(xlim = c(0, max(calibration$Brier_score) * 1.22), clip = "off") +
  labs(title = "Brier score", x = "Lower is better", y = NULL) +
  theme(plot.margin = margin(5.5, 15, 5.5, 5.5))

p_cal_dca <- p_dca +
  (p_brier + theme(legend.position = "none")) +
  plot_layout(widths = c(1.45, 0.9))

save_pub(p_cal_dca, "Fig4d_calibration_decision_curve_summary", 178, 82)

message("Exported Fig4 clinical review panels to: ", out_dir)
