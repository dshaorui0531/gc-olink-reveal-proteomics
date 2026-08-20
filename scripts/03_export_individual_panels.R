#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(pROC)
  library(svglite)
  library(ragg)
})

set.seed(20260819)
options(stringsAsFactors = FALSE)

input_dir <- Sys.getenv("GC_OLINK_SUMMARY_DIR", "data/summary")
result_dir <- Sys.getenv("GC_OLINK_RESULTS_DIR", "results")
tab_dir <- file.path(result_dir, "tables")
panel_dir <- file.path(result_dir, "individual_panels_pub")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = file.path(result_dir, "cache"))

palette <- c(
  H = "#4B5563",
  E = "#3182BD",
  L = "#D24B40",
  normal = "#D8D8D8",
  stage1 = "#D8D8D8",
  stage2 = "#9ECAE1",
  stage3 = "#3182BD",
  stage4 = "#D24B40",
  teal = "#33B5A5",
  orange = "#E28E2C",
  dark = "#272727",
  mid = "#767676"
)

group_labels <- c(H = "Healthy", E = "T1-T2", L = "T3-T4")
group_colors <- c(Healthy = palette[["H"]], `T1-T2` = palette[["E"]], `T3-T4` = palette[["L"]])
comparison_labels <- c(
  Early_vs_Healthy = "T1-T2 vs healthy",
  Late_vs_Early = "T3-T4 vs T1-T2",
  Late_vs_Healthy = "T3-T4 vs healthy",
  Monotonic_trend = "Monotonic trend"
)

figure_font <- "Helvetica"

theme_pub <- function(base_size = 7.2, base_family = figure_font) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(family = base_family, colour = "black", lineheight = 0.95),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.4),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.6, family = base_family, lineheight = 0.95, margin = margin(r = 2)),
      legend.key.size = unit(3.0, "mm"),
      legend.key.width = unit(5.0, "mm"),
      legend.spacing.x = unit(3.0, "mm"),
      legend.spacing.y = unit(1.2, "mm"),
      legend.box.spacing = unit(1.8, "mm"),
      strip.text = element_text(size = base_size - 0.1, face = "bold"),
      plot.title = element_text(size = base_size + 0.7, face = "bold", margin = margin(b = 5)),
      plot.subtitle = element_text(size = base_size - 0.3, margin = margin(b = 5)),
      plot.margin = margin(10, 12, 10, 12),
      panel.grid = element_blank()
    )
}
theme_set(theme_pub())

wrap_text <- function(x, width = 34) {
  vapply(as.character(x), function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

save_panel <- function(plot, name, width_mm = 90, height_mm = 70, dpi = 600) {
  stem <- file.path(panel_dir, name)
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  message("Saving ", name)
  pdf_device <- function(filename, ...) grDevices::pdf(file = filename, ..., family = figure_font, useDingbats = FALSE)
  svg_device <- function(...) svglite::svglite(..., fix_text_size = FALSE, system_fonts = list(sans = figure_font, serif = figure_font, mono = figure_font))
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = w, height = h, device = pdf_device, bg = "white", limitsize = FALSE)
  svg_ok <- try(
    ggplot2::ggsave(paste0(stem, ".svg"), plot, width = w, height = h, device = svg_device, bg = "white", limitsize = FALSE),
    silent = TRUE
  )
  if (inherits(svg_ok, "try-error")) {
    warning("SVG export failed for ", name, "; PDF and TIFF were still exported.")
  }
  ggplot2::ggsave(paste0(stem, ".tiff"), plot, width = w, height = h, dpi = dpi, device = ragg::agg_tiff, bg = "white", limitsize = FALSE)
}

auc_ci_lookup <- function(summary_df, keys) {
  hit <- summary_df$AUC_CI[summary_df$Model %in% keys]
  if (length(hit) == 0 || is.na(hit[1])) return("AUC not available")
  hit[1]
}

read_olink_matrix <- function() {
  npx_path <- file.path(input_dir, "04.Diff_analysis/COND1/all_NPX.xlsx")
  npx <- readxl::read_excel(npx_path, sheet = 1, .name_repair = "unique") |> as.data.frame()
  names(npx)[1] <- "SampleID"
  rownames(npx) <- npx$SampleID
  protein_cols <- setdiff(names(npx), c("SampleID", "Group"))
  expr <- as.matrix(npx[, protein_cols])
  storage.mode(expr) <- "numeric"
  meta <- data.frame(
    SampleID = rownames(expr),
    Group = factor(npx$Group, levels = c("H", "E", "L")),
    GroupShort = factor(group_labels[as.character(npx$Group)], levels = unname(group_labels)),
    stringsAsFactors = FALSE
  )
  list(expr = expr, meta = meta)
}

make_volcano <- function(tt, title, file_tag) {
  df <- tt |>
    mutate(
      neg_log10_FDR = -log10(adj.P.Val),
      regulation = case_when(
        adj.P.Val < 0.05 & logFC > 0 ~ "Up",
        adj.P.Val < 0.05 & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      ),
      label = ifelse(adj.P.Val < 0.05 & rank(adj.P.Val, ties.method = "first") <= 10, Protein, NA)
    )
  label_df <- df |> filter(!is.na(label))
  p <- ggplot(df, aes(logFC, neg_log10_FDR)) +
    geom_point(aes(colour = regulation), size = 1.0, alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.25, colour = palette["mid"]) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.25, colour = palette["mid"]) +
    ggrepel::geom_label_repel(
      data = label_df,
      aes(label = label),
      size = 2.05,
      min.segment.length = 0,
      max.overlaps = Inf,
      box.padding = 0.25,
      point.padding = 0.15,
      label.padding = unit(0.08, "lines"),
      label.r = unit(0.04, "lines"),
      label.size = 0.12,
      fill = ggplot2::alpha("white", 0.86),
      show.legend = FALSE
    ) +
    scale_colour_manual(values = c(Up = "#D24B40", Down = "#3182BD", NS = "#D8D8D8")) +
    coord_cartesian(clip = "off") +
    labs(title = title, x = "log2 fold change", y = "-log10 FDR", colour = NULL) +
    theme(legend.position = "bottom")
  save_panel(p, file_tag, 94, 76)
}

make_heatmap <- function(expr, meta, proteins, title, file_tag, width_mm = 120, height_mm = 95) {
  proteins <- intersect(proteins, colnames(expr))
  y_text_size <- dplyr::case_when(
    length(proteins) > 32 ~ 4.4,
    length(proteins) > 24 ~ 4.9,
    TRUE ~ 5.5
  )
  z <- scale(expr[, proteins, drop = FALSE])
  z[z > 2.5] <- 2.5
  z[z < -2.5] <- -2.5
  ord <- order(meta$Group)
  df <- as.data.frame(z[ord, , drop = FALSE])
  df$SampleID <- rownames(z)[ord]
  long <- pivot_longer(df, cols = all_of(proteins), names_to = "Protein", values_to = "Z")
  long$Protein <- factor(long$Protein, levels = rev(proteins))
  long$SampleID <- factor(long$SampleID, levels = unique(df$SampleID))
  p <- ggplot(long, aes(SampleID, Protein, fill = Z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish) +
    labs(title = title, x = NULL, y = NULL, fill = "Z") +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, barheight = unit(24, "mm"), barwidth = unit(3.0, "mm"))) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = y_text_size, lineheight = 0.88),
      legend.position = "right",
      legend.title = element_text(size = 6.0),
      legend.text = element_text(size = 5.6),
      legend.box.margin = margin(l = 2, r = 4, unit = "mm"),
      plot.margin = margin(t = 2, r = 8, b = 2, l = 2, unit = "mm")
    )
  save_panel(p, file_tag, width_mm, height_mm)
}

dat <- read_olink_matrix()
expr <- dat$expr
meta <- dat$meta

limma_eh <- read.csv(file.path(tab_dir, "limma_E_vs_H.csv"))
limma_le <- read.csv(file.path(tab_dir, "limma_L_vs_E.csv"))
limma_lh <- read.csv(file.path(tab_dir, "limma_L_vs_H.csv"))
trend <- read.csv(file.path(tab_dir, "stage_trend_analysis.csv"))
go_res <- read.csv(file.path(tab_dir, "go_bp_enrichment.csv"))
model_summary <- read.csv(file.path(tab_dir, "diagnostic_model_summary.csv"))
early_pred <- read.csv(file.path(tab_dir, "model_early_vs_healthy_predictions.csv"))
late_pred <- read.csv(file.path(tab_dir, "model_late_vs_early_predictions.csv"))
early_freq <- read.csv(file.path(tab_dir, "model_early_vs_healthy_loocv_feature_frequency.csv"))
late_freq <- read.csv(file.path(tab_dir, "model_late_vs_early_loocv_feature_frequency.csv"))

# Fig. 1 individual panels.
p_count <- meta |>
  count(GroupShort) |>
  ggplot(aes(GroupShort, n, fill = GroupShort)) +
  geom_col(width = 0.62, colour = "black", linewidth = 0.2) +
  geom_text(aes(y = n + 1.3, label = n), vjust = 0, size = 2.7) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(limits = c(0, 37), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Cohort design", x = NULL, y = "Samples") +
  theme(legend.position = "none")
save_panel(p_count, "Fig1a_cohort_design", 82, 72)

pca <- prcomp(expr, center = TRUE, scale. = TRUE)
pca_df <- data.frame(SampleID = rownames(expr), PC1 = pca$x[, 1], PC2 = pca$x[, 2], Group = meta$GroupShort)
var_exp <- (pca$sdev^2 / sum(pca$sdev^2))[1:2] * 100
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = Group)) +
  geom_point(size = 2.2, alpha = 0.9) +
  stat_ellipse(linewidth = 0.35, level = 0.68) +
  scale_colour_manual(values = group_colors) +
  labs(title = "Global protein profile", x = sprintf("PC1 (%.1f%%)", var_exp[1]), y = sprintf("PC2 (%.1f%%)", var_exp[2]), colour = NULL) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7.0, family = figure_font, margin = margin(r = 5)),
    plot.margin = margin(10, 12, 14, 12)
  ) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 2.4)))
save_panel(p_pca, "Fig1b_pca_global_profile", 90, 76)

density_df <- as.data.frame(expr) |>
  mutate(SampleID = rownames(expr), Group = meta$GroupShort) |>
  pivot_longer(cols = -c(SampleID, Group), names_to = "Protein", values_to = "NPX")
p_density <- ggplot(density_df, aes(NPX, colour = Group)) +
  geom_density(linewidth = 0.45, key_glyph = "path") +
  scale_colour_manual(values = group_colors) +
  labs(title = "NPX distribution", x = "NPX", y = "Density", colour = NULL) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7.0, family = figure_font, margin = margin(r = 5)),
    plot.margin = margin(10, 12, 14, 12)
  ) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = 0.85)))
save_panel(p_density, "Fig1c_npx_distribution", 90, 74)

cor_mat <- cor(t(expr), method = "pearson")
cor_long <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
names(cor_long) <- c("Sample1", "Sample2", "Correlation")
cor_long <- cor_long |>
  filter(as.character(Sample1) < as.character(Sample2)) |>
  left_join(meta |> dplyr::select(SampleID, Group1 = GroupShort), by = c("Sample1" = "SampleID")) |>
  left_join(meta |> dplyr::select(SampleID, Group2 = GroupShort), by = c("Sample2" = "SampleID")) |>
  mutate(Pair = ifelse(Group1 == Group2, paste0("Within ", Group1), "Between groups"))
p_cor <- ggplot(cor_long, aes(Pair, Correlation)) +
  geom_boxplot(outlier.size = 0.35, width = 0.6, fill = "#D8D8D8", colour = "black", linewidth = 0.25) +
  labs(title = "Sample-level correlation", x = NULL, y = "Pearson r") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_panel(p_cor, "Fig1d_sample_correlation", 95, 74)

# Fig. 2 individual panels.
make_volcano(limma_eh, "T1-T2 GC vs healthy", "Fig2a_volcano_E_vs_H")
make_volcano(limma_le, "T3-T4 vs T1-T2 GC", "Fig2b_volcano_L_vs_E")
make_volcano(limma_lh, "T3-T4 GC vs healthy", "Fig2c_volcano_L_vs_H")

sig_count <- bind_rows(
  limma_eh |> mutate(Comparison = "E vs H"),
  limma_le |> mutate(Comparison = "L vs E"),
  limma_lh |> mutate(Comparison = "L vs H")
) |>
  mutate(Comparison = recode(Comparison, `E vs H` = "T1-T2 vs H", `L vs E` = "T3-T4 vs T1-T2", `L vs H` = "T3-T4 vs H")) |>
  mutate(Direction = case_when(adj.P.Val < 0.05 & logFC > 0 ~ "Up",
                               adj.P.Val < 0.05 & logFC < 0 ~ "Down",
                               TRUE ~ "NS")) |>
  filter(Direction != "NS") |>
  count(Comparison, Direction)
p_sigcount <- ggplot(sig_count, aes(Comparison, n, fill = Direction)) +
  geom_col(width = 0.65, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c(Up = "#D24B40", Down = "#3182BD")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "FDR < 0.05 proteins", x = NULL, y = "Count", fill = NULL) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6.8),
    legend.key.width = unit(5.2, "mm"),
    legend.spacing.x = unit(1.2, "mm")
  )
save_panel(p_sigcount, "Fig2d_differential_count", 88, 68)

early_sig <- limma_eh |> filter(adj.P.Val < 0.05) |> arrange(adj.P.Val)
late_sig <- limma_le |> filter(adj.P.Val < 0.05) |> arrange(adj.P.Val)
hl_sig <- limma_lh |> filter(adj.P.Val < 0.05) |> arrange(adj.P.Val)
top_heat_proteins <- unique(c(head(early_sig$Protein, 8), head(late_sig$Protein, 8), head(hl_sig$Protein, 10)))
top_heat_proteins <- top_heat_proteins[seq_len(min(24, length(top_heat_proteins)))]
make_heatmap(expr, meta, top_heat_proteins, "Top differential proteins", "Fig2e_top_differential_heatmap", 128, 110)

# Fig. 3 individual panels.
p_trend_vol_df <- trend |>
  mutate(Signature = case_when(monotonic & trend_FDR < 0.05 & slope > 0 ~ "Increasing",
                               monotonic & trend_FDR < 0.05 & slope < 0 ~ "Decreasing",
                               TRUE ~ "Other"),
         label = ifelse(Signature != "Other" & rank(trend_FDR, ties.method = "first") <= 12, Protein, NA))
p_trend_vol_label <- p_trend_vol_df |> filter(!is.na(label))
p_trend_vol <- p_trend_vol_df |>
  ggplot(aes(slope, -log10(trend_FDR), colour = Signature)) +
  geom_point(size = 0.95, alpha = 0.85) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.25, colour = palette["mid"]) +
  ggrepel::geom_label_repel(
    data = p_trend_vol_label,
    aes(label = label),
    size = 2.05,
    min.segment.length = 0,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
    label.padding = unit(0.08, "lines"),
    label.r = unit(0.04, "lines"),
    label.size = 0.12,
    fill = ggplot2::alpha("white", 0.86),
    colour = "black",
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c(Increasing = "#D24B40", Decreasing = "#3182BD", Other = "#D8D8D8")) +
  coord_cartesian(clip = "off") +
  labs(title = "Ordered T-category trend", x = "Slope per group step", y = "-log10 trend FDR", colour = NULL) +
  theme(legend.position = "bottom")
save_panel(p_trend_vol, "Fig3a_trend_volcano", 96, 78)

trend_sig <- trend |> filter(monotonic, trend_FDR < 0.05) |> arrange(trend_FDR)
trend_plot_proteins <- head(trend_sig$Protein, 12)
trend_line_df <- as.data.frame(expr[, trend_plot_proteins, drop = FALSE]) |>
  mutate(SampleID = rownames(expr), Group = meta$GroupShort) |>
  pivot_longer(cols = all_of(trend_plot_proteins), names_to = "Protein", values_to = "NPX") |>
  mutate(Protein = factor(Protein, levels = trend_plot_proteins))
p_trend_line <- ggplot(trend_line_df, aes(Group, NPX, group = Protein, colour = Protein)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.45, alpha = 0.9) +
  stat_summary(fun = mean, geom = "point", size = 1.5) +
  labs(title = "Top monotonic proteins", x = NULL, y = "Mean NPX") +
  theme(legend.position = "right", legend.key.height = unit(3.8, "mm"))
save_panel(p_trend_line, "Fig3b_top_monotonic_lineplot", 105, 82)

trend_top_heat <- head(trend_sig$Protein, 30)
make_heatmap(expr, meta, trend_top_heat, "Monotonic T-category protein signature", "Fig3c_monotonic_heatmap", 150, 132)

# Fig. 4 individual panels.
roc_early <- pROC::roc(early_pred$Outcome, early_pred$CV_Probability, quiet = TRUE, levels = c(0, 1), direction = "<")
roc_late <- pROC::roc(late_pred$Outcome, late_pred$CV_Probability, quiet = TRUE, levels = c(0, 1), direction = "<")
roc_df <- bind_rows(
  data.frame(FPR = 1 - roc_early$specificities, Sensitivity = roc_early$sensitivities, Model = "T1-T2 GC vs healthy"),
  data.frame(FPR = 1 - roc_late$specificities, Sensitivity = roc_late$sensitivities, Model = "T3-T4 vs T1-T2 GC")
)
p_roc <- ggplot(roc_df, aes(FPR, Sensitivity, colour = Model)) +
  geom_line(linewidth = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.25, colour = palette["mid"]) +
  coord_equal() +
  scale_colour_manual(
    values = c("T1-T2 GC vs healthy" = "#3182BD", "T3-T4 vs T1-T2 GC" = "#D24B40"),
    labels = c(
      paste0("T1-T2 model: ", auc_ci_lookup(model_summary, c("Early GC vs healthy", "T1-T2 GC vs healthy"))),
      paste0("Invasion-depth model: ", auc_ci_lookup(model_summary, c("Late GC vs early GC", "T3-T4 vs T1-T2 GC")))
    )
  ) +
  labs(title = "Nested LOOCV ROC",
       x = "1 - specificity", y = "Sensitivity") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(7, "mm"),
    legend.text = element_text(size = 6.2)
  ) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))
save_panel(p_roc, "Fig4a_nested_loocv_roc", 96, 96)

plot_pred <- function(df, title, name) {
  df$GroupLabel <- factor(group_labels[as.character(df$Group)], levels = unname(group_labels))
  p <- ggplot(df, aes(GroupLabel, CV_Probability, fill = GroupLabel)) +
    geom_boxplot(width = 0.58, outlier.size = 0.6, colour = "black", linewidth = 0.25) +
    geom_jitter(width = 0.08, size = 0.8, alpha = 0.75) +
    scale_fill_manual(values = group_colors) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.06))) +
    labs(title = title, x = NULL, y = "Probability of positive class") +
    theme(legend.position = "none")
  save_panel(p, name, 78, 72)
}
plot_pred(early_pred, "T1-T2 GC vs healthy", "Fig4b_predicted_probability_E_vs_H")
plot_pred(late_pred, "T3-T4 vs T1-T2 GC", "Fig4c_predicted_probability_L_vs_E")

plot_freq <- function(df, title, name, fill) {
  df <- df |> arrange(desc(Frequency), Protein) |> slice_head(n = min(12, nrow(df)))
  p <- ggplot(df, aes(reorder(Protein, Frequency), Frequency)) +
    geom_col(width = 0.7, fill = fill, colour = "black", linewidth = 0.15) +
    coord_flip() +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.08))) +
    labs(title = title, x = NULL, y = "LOOCV selection frequency")
  save_panel(p, name, 92, 86)
}
plot_freq(early_freq, "Feature stability: T1-T2 model", "Fig4d_feature_stability_early", "#3182BD")
plot_freq(late_freq, "Feature stability: invasion-depth model", "Fig4e_feature_stability_late", "#D24B40")

# Fig. 5 individual panels by enrichment gene set.
if (nrow(go_res) > 0) {
  for (gs in unique(go_res$GeneSet)) {
    gdf <- go_res |>
      filter(GeneSet == gs) |>
      arrange(p.adjust) |>
      slice_head(n = min(10, nrow(filter(go_res, GeneSet == gs)))) |>
      mutate(DescriptionWrapped = factor(wrap_text(Description, 36), levels = rev(wrap_text(Description, 36))))
    p <- ggplot(gdf, aes(-log10(p.adjust), DescriptionWrapped)) +
      geom_point(colour = "#33B5A5", size = 2.2, alpha = 0.9) +
      scale_x_continuous(expand = expansion(mult = c(0.03, 0.16))) +
      labs(title = paste("GO BP:", ifelse(gs %in% names(comparison_labels), comparison_labels[[gs]], gsub("_", " ", gs))), x = "-log10(FDR)", y = NULL) +
      theme(
        axis.text.y = element_text(size = 6.2, family = figure_font, lineheight = 0.9, margin = margin(r = 3)),
        axis.text.x = element_text(family = figure_font),
        axis.title.x = element_text(family = figure_font),
        plot.title = element_text(family = figure_font),
        legend.position = "none"
      )
    p <- p + theme(plot.margin = margin(t = 3, r = 8, b = 4, l = 14, unit = "mm"))
    save_panel(p, paste0("Fig5_", gs, "_GO_BP"), 150, 98)
  }
}

# Fig. 6 individual panels.
tcga_expr <- read.csv(file.path(tab_dir, "tcga_stad_candidate_expression_long.csv"))
tcga_stage <- read.csv(file.path(tab_dir, "tcga_stad_stage_association.csv"))
tcga_stage_expr <- tcga_expr |> filter(Tissue == "Tumor")

stage_tbl <- read.csv(file.path(tab_dir, "tcga_stad_stage_association.csv"))
surv_tbl <- read.csv(file.path(tab_dir, "tcga_stad_candidate_survival.csv"))
sig_scores <- read.csv(file.path(tab_dir, "tcga_stad_signature_scores.csv"))

tn_genes <- c("CCN1", "DSG3", "ITGA11", "KDR", "MUC16", "NPY", "NXPH3", "RET", "ROBO2", "TNFRSF10B")
for (g in intersect(tn_genes, unique(tcga_expr$Symbol))) {
  gdf <- tcga_expr |>
    filter(Symbol == g, Tissue %in% c("Adjacent normal", "Tumor")) |>
    mutate(TissueShort = factor(ifelse(Tissue == "Adjacent normal", "Normal", "Tumor"), levels = c("Normal", "Tumor")))
  p <- ggplot(gdf, aes(TissueShort, TPM_log2, fill = TissueShort)) +
    geom_boxplot(width = 0.58, outlier.size = 0.35, colour = "black", linewidth = 0.25) +
    scale_fill_manual(values = c(Normal = "#D8D8D8", Tumor = "#D24B40")) +
    labs(title = paste0("TCGA-STAD: ", g), x = NULL, y = "STAR TPM (log2)") +
    theme(legend.position = "none")
  save_panel(p, paste0("Fig6a_tumor_normal_", g), 68, 72)
}

# Reconstruct stage labels from the saved long table is not possible alone; use the
# full source table created by the validation script when available.
stage_source_path <- file.path(tab_dir, "tcga_stad_stage_expression_long.csv")
if (file.exists(stage_source_path)) {
  stage_expr <- read.csv(stage_source_path)
  stage_genes <- stage_tbl |> arrange(spearman_p) |> slice_head(n = min(8, nrow(stage_tbl))) |> pull(Symbol)
  for (g in intersect(stage_genes, unique(stage_expr$Symbol))) {
    gdf <- stage_expr |>
      filter(Symbol == g) |>
      mutate(Stage = factor(Stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
             StageShort = factor(gsub("Stage ", "", Stage), levels = c("I", "II", "III", "IV")))
    p <- ggplot(gdf, aes(StageShort, TPM_log2, fill = StageShort)) +
      geom_boxplot(width = 0.58, outlier.size = 0.35, colour = "black", linewidth = 0.25) +
      scale_fill_manual(values = c("I" = "#D8D8D8", "II" = "#9ECAE1", "III" = "#3182BD", "IV" = "#D24B40")) +
      labs(title = paste0("Stage association: ", g), x = "Pathological stage", y = "STAR TPM (log2)") +
      theme(legend.position = "none")
    save_panel(p, paste0("Fig6b_stage_", g), 70, 72)
  }
}

forest_df <- surv_tbl |>
  filter(!is.na(cox_p)) |>
  arrange(cox_p) |>
  slice_head(n = min(14, nrow(surv_tbl))) |>
  mutate(Gene = factor(Gene, levels = rev(Gene)))
p_forest <- ggplot(forest_df, aes(HR_high_vs_low, Gene)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#767676", linewidth = 0.25) +
  geom_errorbar(aes(xmin = HR_lower95, xmax = HR_upper95), orientation = "y", width = 0.18, linewidth = 0.35) +
  geom_point(size = 1.6, colour = "#272727") +
  scale_x_log10() +
  labs(title = "Overall survival association", x = "Hazard ratio, high vs low expression", y = NULL)
save_panel(p_forest, "Fig6c_survival_forest", 108, 90)

for (sig in unique(sig_scores$Signature)) {
  sdf <- sig_scores |>
    filter(Signature == sig) |>
    mutate(Group = ifelse(Score > median(Score, na.rm = TRUE), "High", "Low"),
           Group = factor(Group, levels = c("High", "Low")))
  p <- ggplot(sdf, aes(Group, Score, fill = Group)) +
    geom_boxplot(width = 0.58, outlier.size = 0.45, colour = "black", linewidth = 0.25) +
    scale_fill_manual(values = c(High = "#D24B40", Low = "#3182BD")) +
    labs(
      title = recode(
        sig,
        `Early detection signature` = "T1-T2 detection signature",
        `Early-detection signature` = "T1-T2 detection signature",
        `Early_detection_signature` = "T1-T2 detection signature",
        `Progression signature` = "Invasion-depth signature",
        `Progression_signature` = "Invasion-depth signature"
      ),
      x = NULL,
      y = "Coefficient-weighted z score"
    ) +
    theme(legend.position = "none")
  save_panel(p, paste0("Fig6d_signature_score_", gsub("[^A-Za-z0-9]+", "_", sig)), 82, 72)
}

message("Individual panel export complete: ", panel_dir)
