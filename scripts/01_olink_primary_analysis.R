#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(limma)
  library(pROC)
  library(glmnet)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(svglite)
  library(ragg)
})

set.seed(20260819)
options(stringsAsFactors = FALSE)

input_dir <- Sys.getenv("GC_OLINK_SUMMARY_DIR", "data/summary")
out_dir <- Sys.getenv("GC_OLINK_RESULTS_DIR", "results")
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
report_dir <- file.path(out_dir, "reports")
cache_dir <- file.path(out_dir, "cache")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = cache_dir)

palette <- c(
  H = "#4B5563",
  E = "#3182BD",
  L = "#D24B40",
  neutral_dark = "#272727",
  neutral_mid = "#767676",
  neutral_light = "#D8D8D8",
  signal_teal = "#33B5A5",
  accent_orange = "#E28E2C"
)

figure_font <- "Helvetica"
group_labels <- c(H = "Healthy", E = "T1-T2", L = "T3-T4")
group_colors <- c(Healthy = palette[["H"]], `T1-T2` = palette[["E"]], `T3-T4` = palette[["L"]])
model_colors <- c(`T1-T2 GC vs healthy` = "#3182BD", `T3-T4 vs T1-T2 GC` = "#D24B40")

theme_pub <- function(base_size = 7, base_family = figure_font) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(family = base_family, colour = "black", lineheight = 0.95),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.6, family = base_family, lineheight = 0.95, margin = margin(r = 2)),
      legend.key.size = unit(3.0, "mm"),
      legend.spacing.x = unit(1.0, "mm"),
      legend.spacing.y = unit(1.0, "mm"),
      strip.text = element_text(size = base_size - 0.2, face = "bold"),
      plot.title = element_text(size = base_size + 0.6, face = "bold"),
      plot.margin = margin(t = 3, r = 6, b = 3, l = 3, unit = "mm"),
      panel.grid = element_blank()
    )
}
theme_set(theme_pub())

save_pub <- function(plot, stem, width_mm = 183, height_mm = 130, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  message("Saving figure: ", basename(stem))
  message("  - PDF")
  pdf_device <- function(filename, ...) grDevices::pdf(file = filename, ..., family = figure_font, useDingbats = FALSE)
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = w, height = h, device = pdf_device)
  message("  - TIFF")
  ggplot2::ggsave(paste0(stem, ".tiff"), plot, width = w, height = h, dpi = dpi, device = ragg::agg_tiff)
}

safe_write <- function(x, file) {
  write.csv(x, file, row.names = FALSE, fileEncoding = "UTF-8")
}

auc_ci_text <- function(roc_obj) {
  ci <- as.numeric(pROC::ci.auc(roc_obj))
  sprintf("AUC %.3f (95%% CI %.3f-%.3f)", as.numeric(pROC::auc(roc_obj)), ci[1], ci[3])
}

read_input_data <- function() {
  npx_path <- file.path(input_dir, "04.Diff_analysis/COND1/all_NPX.xlsx")
  anno_path <- file.path(input_dir, "03.Profile/anno_profile.xlsx")
  sample_path <- file.path(input_dir, "01.Raw/Sample_Info.txt")

  npx <- readxl::read_excel(npx_path, sheet = 1, .name_repair = "unique") |> as.data.frame()
  names(npx)[1] <- "SampleID"
  rownames(npx) <- npx$SampleID
  protein_cols <- setdiff(names(npx), c("SampleID", "Group"))
  expr <- as.matrix(npx[, protein_cols])
  storage.mode(expr) <- "numeric"

  group_from_matrix <- if ("Group" %in% names(npx)) npx$Group else substr(rownames(expr), 1, 1)
  meta <- data.frame(
    SampleID = rownames(expr),
    Group = factor(group_from_matrix, levels = c("H", "E", "L")),
    stringsAsFactors = FALSE
  )
  meta$GroupLabel <- factor(
    group_labels[as.character(meta$Group)],
    levels = unname(group_labels)
  )

  raw_sample_info <- read.delim(sample_path, check.names = FALSE)
  anno <- readxl::read_excel(anno_path, sheet = 1, .name_repair = "unique") |> as.data.frame()
  anno <- anno |> distinct(Assay, .keep_all = TRUE)
  list(expr = expr, meta = meta, sample_info = raw_sample_info, anno = anno)
}

impute_median <- function(mat) {
  out <- mat
  for (j in seq_len(ncol(out))) {
    miss <- is.na(out[, j])
    if (any(miss)) out[miss, j] <- median(out[, j], na.rm = TRUE)
  }
  out
}

run_limma <- function(expr, meta) {
  design <- model.matrix(~ 0 + Group, data = meta)
  colnames(design) <- levels(meta$Group)
  fit <- lmFit(t(expr), design)
  cont <- makeContrasts(
    E_vs_H = E - H,
    L_vs_E = L - E,
    L_vs_H = L - H,
    levels = design
  )
  fit2 <- eBayes(contrasts.fit(fit, cont))
  lapply(colnames(cont), function(cn) {
    tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "P")
    tt$Protein <- rownames(tt)
    tt <- tt[, c("Protein", setdiff(names(tt), "Protein"))]
    tt
  }) |> setNames(colnames(cont))
}

run_trend <- function(expr, meta) {
  gnum <- dplyr::recode(as.character(meta$Group), H = 0, E = 1, L = 2) |> as.numeric()
  res <- lapply(colnames(expr), function(p) {
    x <- expr[, p]
    means <- tapply(x, meta$Group, mean, na.rm = TRUE)
    fit <- summary(lm(x ~ gnum))
    an <- summary(aov(x ~ meta$Group))[[1]]
    data.frame(
      Protein = p,
      slope = unname(coef(fit)[2, "Estimate"]),
      trend_p = unname(coef(fit)[2, "Pr(>|t|)"]),
      anova_p = an[1, "Pr(>F)"],
      mean_H = unname(means["H"]),
      mean_E = unname(means["E"]),
      mean_L = unname(means["L"]),
      monotonic = (means["H"] < means["E"] && means["E"] < means["L"]) ||
        (means["H"] > means["E"] && means["E"] > means["L"]),
      direction = ifelse(means["H"] < means["E"] && means["E"] < means["L"], "increasing",
                         ifelse(means["H"] > means["E"] && means["E"] > means["L"], "decreasing", "non-monotonic"))
    )
  }) |> bind_rows()
  res$trend_FDR <- p.adjust(res$trend_p, method = "BH")
  res$anova_FDR <- p.adjust(res$anova_p, method = "BH")
  res |> arrange(trend_FDR, anova_FDR)
}

run_go_enrichment <- function(gene_sets, universe_symbols) {
  universe_map <- suppressMessages(
    bitr(unique(universe_symbols), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )
  out <- list()
  for (nm in names(gene_sets)) {
    gm <- suppressMessages(
      bitr(unique(gene_sets[[nm]]), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    )
    if (nrow(gm) < 5) next
    eg <- suppressMessages(enrichGO(
      gene = unique(gm$ENTREZID),
      universe = unique(universe_map$ENTREZID),
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    ))
    if (!is.null(eg) && nrow(as.data.frame(eg)) > 0) {
      df <- as.data.frame(eg)
      df$GeneSet <- nm
      out[[nm]] <- df
    }
  }
  if (length(out) == 0) return(data.frame())
  bind_rows(out) |> arrange(GeneSet, p.adjust)
}

rank_by_training_ttest <- function(x_train, y_train, top_n = 80) {
  pvals <- apply(x_train, 2, function(v) {
    out <- try(t.test(v ~ y_train)$p.value, silent = TRUE)
    if (inherits(out, "try-error") || is.na(out)) 1 else out
  })
  names(sort(pvals))[seq_len(min(top_n, length(pvals)))]
}

nested_glmnet_binary <- function(expr, meta, positive, negative, label, top_n = 5) {
  keep <- meta$Group %in% c(negative, positive)
  x <- expr[keep, , drop = FALSE]
  m <- droplevels(meta[keep, ])
  y <- ifelse(m$Group == positive, 1, 0)
  pred <- rep(NA_real_, length(y))
  selected <- list()
  fold_lambda <- rep(NA_real_, length(y))

  for (i in seq_along(y)) {
    train <- setdiff(seq_along(y), i)
    feats <- rank_by_training_ttest(x[train, , drop = FALSE], y[train], top_n = top_n)
    inner_nfolds <- min(5, min(table(y[train])))
    cvfit <- cv.glmnet(
      x = x[train, feats, drop = FALSE],
      y = y[train],
      family = "binomial",
      alpha = 0,
      nfolds = inner_nfolds,
      type.measure = "deviance",
      standardize = TRUE
    )
    s_use <- "lambda.1se"
    pred[i] <- as.numeric(predict(cvfit, newx = x[i, feats, drop = FALSE], s = s_use, type = "response"))
    fold_lambda[i] <- cvfit[[s_use]]
    cc <- as.matrix(coef(cvfit, s = s_use))
    selected[[i]] <- rownames(cc)[which(cc[, 1] != 0)]
    selected[[i]] <- setdiff(selected[[i]], "(Intercept)")
  }

  roc_obj <- pROC::roc(y, pred, quiet = TRUE, levels = c(0, 1), direction = "<")

  final_feats <- rank_by_training_ttest(x, y, top_n = top_n)
  final_nfolds <- min(10, min(table(y)))
  final_cv <- cv.glmnet(
    x = x[, final_feats, drop = FALSE],
    y = y,
    family = "binomial",
    alpha = 0,
    nfolds = final_nfolds,
    type.measure = "deviance",
    standardize = TRUE
  )
  final_s <- "lambda.1se"
  final_coef <- as.matrix(coef(final_cv, s = final_s))
  coef_table <- data.frame(
    Protein = rownames(final_coef),
    Coefficient = final_coef[, 1],
    row.names = NULL
  ) |>
    filter(Coefficient != 0)
  final_score <- as.numeric(predict(final_cv, newx = x[, final_feats, drop = FALSE], s = final_s, type = "response"))

  freq <- table(unlist(selected))
  freq_table <- data.frame(
    Protein = names(freq),
    FoldCount = as.integer(freq),
    Frequency = as.integer(freq) / length(y),
    row.names = NULL
  ) |>
    arrange(desc(Frequency), Protein)

  list(
    label = label,
    meta = m,
    y = y,
    pred_cv = pred,
    final_score = final_score,
    roc = roc_obj,
    coef = coef_table,
    feature_frequency = freq_table,
    fold_lambda = fold_lambda
  )
}

plot_roc_df <- function(model_result) {
  data.frame(
    Specificity = model_result$roc$specificities,
    Sensitivity = model_result$roc$sensitivities,
    FPR = 1 - model_result$roc$specificities,
    Model = model_result$label
  )
}

make_volcano <- function(tt, title) {
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
  ggplot(df, aes(logFC, neg_log10_FDR)) +
    geom_point(aes(colour = regulation), size = 0.9, alpha = 0.82) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.25, colour = palette["neutral_mid"]) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.25, colour = palette["neutral_mid"]) +
    ggrepel::geom_text_repel(data = label_df, aes(label = label), size = 2, min.segment.length = 0, max.overlaps = Inf) +
    scale_color_manual(
      values = c(Up = "#D24B40", Down = "#3182BD", NS = "#D8D8D8"),
      breaks = c("Down", "NS", "Up"),
      guide = guide_legend(title = NULL, override.aes = list(size = 2), nrow = 1, byrow = TRUE)
    ) +
    labs(title = title, x = "log2 fold change", y = "-log10 FDR") +
    theme(legend.position = "bottom")
}

make_heatmap <- function(expr, meta, proteins, title, zlim = 2.5) {
  proteins <- intersect(proteins, colnames(expr))
  y_text_size <- dplyr::case_when(
    length(proteins) > 32 ~ 4.4,
    length(proteins) > 24 ~ 4.9,
    TRUE ~ 5.5
  )
  z <- scale(expr[, proteins, drop = FALSE])
  z[z > zlim] <- zlim
  z[z < -zlim] <- -zlim
  ord <- order(meta$Group)
  df <- as.data.frame(z[ord, , drop = FALSE])
  df$SampleID <- rownames(z)[ord]
  df$Group <- meta$Group[ord]
  long <- pivot_longer(df, cols = all_of(proteins), names_to = "Protein", values_to = "Z")
  long$Protein <- factor(long$Protein, levels = rev(proteins))
  long$SampleID <- factor(long$SampleID, levels = unique(df$SampleID))
  ggplot(long, aes(SampleID, Protein, fill = Z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0, limits = c(-zlim, zlim), oob = scales::squish) +
    labs(title = title, x = NULL, y = NULL, fill = "Z") +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, barheight = unit(22, "mm"), barwidth = unit(3.0, "mm"))) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = y_text_size, lineheight = 0.88),
      legend.title = element_text(size = 6.0),
      legend.text = element_text(size = 5.6),
      legend.box.margin = margin(l = 2, r = 4, unit = "mm"),
      legend.position = "right"
    )
}

dat <- read_input_data()
expr_raw <- dat$expr
meta <- dat$meta
anno <- dat$anno
expr <- impute_median(expr_raw)

qc_summary <- data.frame(
  Metric = c("Samples", "Proteins", "Missing values", "Maximum missing fraction per protein", "Groups"),
  Value = c(
    nrow(expr),
    ncol(expr),
    sum(is.na(expr_raw)),
    sprintf("%.3f", max(colMeans(is.na(expr_raw)))),
    paste(names(table(meta$Group)), as.integer(table(meta$Group)), sep = "=", collapse = "; ")
  )
)
safe_write(qc_summary, file.path(tab_dir, "qc_summary.csv"))
safe_write(meta, file.path(tab_dir, "sample_metadata.csv"))

protein_qc <- data.frame(
  Protein = colnames(expr_raw),
  MissingFraction = colMeans(is.na(expr_raw)),
  MeanNPX = colMeans(expr, na.rm = TRUE),
  SdNPX = apply(expr, 2, sd, na.rm = TRUE),
  row.names = NULL
) |>
  arrange(desc(SdNPX))
safe_write(protein_qc, file.path(tab_dir, "protein_qc_summary.csv"))

diff_res <- run_limma(expr, meta)
for (nm in names(diff_res)) {
  safe_write(diff_res[[nm]], file.path(tab_dir, paste0("limma_", nm, ".csv")))
}

sig_sets <- lapply(diff_res, function(x) x$Protein[x$adj.P.Val < 0.05])
venn_members <- data.frame(
  Protein = colnames(expr),
  E_vs_H = colnames(expr) %in% sig_sets$E_vs_H,
  L_vs_E = colnames(expr) %in% sig_sets$L_vs_E,
  L_vs_H = colnames(expr) %in% sig_sets$L_vs_H
) |>
  mutate(SignificantComparisons = E_vs_H + L_vs_E + L_vs_H) |>
  arrange(desc(SignificantComparisons), Protein)
safe_write(venn_members, file.path(tab_dir, "differential_intersections.csv"))

trend_res <- run_trend(expr, meta)
safe_write(trend_res, file.path(tab_dir, "stage_trend_analysis.csv"))

trend_sig <- trend_res |>
  filter(monotonic, trend_FDR < 0.05) |>
  arrange(trend_FDR)
early_sig <- diff_res$E_vs_H |> filter(adj.P.Val < 0.05) |> arrange(adj.P.Val)
late_sig <- diff_res$L_vs_E |> filter(adj.P.Val < 0.05) |> arrange(adj.P.Val)
hl_sig <- diff_res$L_vs_H |> filter(adj.P.Val < 0.05) |> arrange(adj.P.Val)

candidate_summary <- bind_rows(
  early_sig |> transmute(Protein, EvidenceSet = "T1-T2 GC vs healthy", logFC, FDR = adj.P.Val),
  late_sig |> transmute(Protein, EvidenceSet = "T3-T4 vs T1-T2 GC", logFC, FDR = adj.P.Val),
  trend_sig |> transmute(Protein, EvidenceSet = "Monotonic T-category trend", logFC = slope, FDR = trend_FDR)
) |>
  left_join(anno |> dplyr::select(Assay, OlinkID, UniProt, description), by = c("Protein" = "Assay"))
safe_write(candidate_summary, file.path(tab_dir, "candidate_protein_summary.csv"))

gene_sets <- list(
  Early_vs_Healthy = early_sig$Protein,
  Late_vs_Early = late_sig$Protein,
  Late_vs_Healthy = hl_sig$Protein,
  Monotonic_trend = trend_sig$Protein
)
go_res <- run_go_enrichment(gene_sets, colnames(expr))
safe_write(go_res, file.path(tab_dir, "go_bp_enrichment.csv"))

model_early <- nested_glmnet_binary(expr, meta, positive = "E", negative = "H", label = "T1-T2 GC vs healthy")
model_late <- nested_glmnet_binary(expr, meta, positive = "L", negative = "E", label = "T3-T4 vs T1-T2 GC")

safe_write(model_early$coef, file.path(tab_dir, "model_early_vs_healthy_final_coefficients.csv"))
safe_write(model_early$feature_frequency, file.path(tab_dir, "model_early_vs_healthy_loocv_feature_frequency.csv"))
safe_write(data.frame(SampleID = model_early$meta$SampleID, Group = model_early$meta$Group, Outcome = model_early$y, CV_Probability = model_early$pred_cv, FinalModelScore = model_early$final_score),
           file.path(tab_dir, "model_early_vs_healthy_predictions.csv"))

safe_write(model_late$coef, file.path(tab_dir, "model_late_vs_early_final_coefficients.csv"))
safe_write(model_late$feature_frequency, file.path(tab_dir, "model_late_vs_early_loocv_feature_frequency.csv"))
safe_write(data.frame(SampleID = model_late$meta$SampleID, Group = model_late$meta$Group, Outcome = model_late$y, CV_Probability = model_late$pred_cv, FinalModelScore = model_late$final_score),
           file.path(tab_dir, "model_late_vs_early_predictions.csv"))

model_summary <- data.frame(
  Model = c(model_early$label, model_late$label),
  Nested_LOOCV_AUC = c(as.numeric(auc(model_early$roc)), as.numeric(auc(model_late$roc))),
  AUC_CI = c(auc_ci_text(model_early$roc), auc_ci_text(model_late$roc)),
  FinalModelFeatures = c(
    paste(setdiff(model_early$coef$Protein, "(Intercept)"), collapse = ";"),
    paste(setdiff(model_late$coef$Protein, "(Intercept)"), collapse = ";")
  )
)
safe_write(model_summary, file.path(tab_dir, "diagnostic_model_summary.csv"))

# Figure 1: cohort and global NPX structure.
p_count <- meta |>
  count(GroupLabel) |>
  ggplot(aes(GroupLabel, n, fill = GroupLabel)) +
  geom_col(width = 0.65, colour = "black", linewidth = 0.2) +
  geom_text(aes(y = n + 1.2, label = n), vjust = 0, size = 2.4) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(limits = c(0, 36), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Cohort design", x = NULL, y = "Samples") +
  theme(legend.position = "none")

pca <- prcomp(expr, center = TRUE, scale. = TRUE)
pca_df <- data.frame(
  SampleID = rownames(expr),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = meta$GroupLabel
)
var_exp <- (pca$sdev^2 / sum(pca$sdev^2))[1:2] * 100
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = Group)) +
  geom_point(size = 2.2, alpha = 0.9) +
  stat_ellipse(linewidth = 0.35, level = 0.68) +
  scale_colour_manual(values = group_colors) +
  labs(title = "Global protein profile", x = sprintf("PC1 (%.1f%%)", var_exp[1]), y = sprintf("PC2 (%.1f%%)", var_exp[2])) +
  theme(legend.position = "bottom")

density_df <- as.data.frame(expr) |>
  mutate(SampleID = rownames(expr), Group = meta$GroupLabel) |>
  pivot_longer(cols = -c(SampleID, Group), names_to = "Protein", values_to = "NPX")
p_density <- ggplot(density_df, aes(NPX, colour = Group)) +
  geom_density(linewidth = 0.45) +
  scale_colour_manual(values = group_colors) +
  labs(title = "NPX distribution", x = "NPX", y = "Density") +
  theme(legend.position = "bottom")

cor_mat <- cor(t(expr), method = "pearson")
cor_long <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
names(cor_long) <- c("Sample1", "Sample2", "Correlation")
cor_long <- cor_long |>
  filter(as.character(Sample1) < as.character(Sample2)) |>
  left_join(meta |> dplyr::select(SampleID, Group1 = GroupLabel), by = c("Sample1" = "SampleID")) |>
  left_join(meta |> dplyr::select(SampleID, Group2 = GroupLabel), by = c("Sample2" = "SampleID")) |>
  mutate(Pair = ifelse(Group1 == Group2, paste0("Within ", Group1), "Between groups"))
p_cor <- ggplot(cor_long, aes(Pair, Correlation)) +
  geom_boxplot(outlier.size = 0.35, width = 0.6, fill = "#D8D8D8", colour = "black", linewidth = 0.25) +
  labs(title = "Sample-level correlation", x = NULL, y = "Pearson r") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

fig1 <- ((p_count | p_pca) / (p_density | p_cor)) +
  plot_annotation(tag_levels = "a")
save_pub(fig1, file.path(fig_dir, "Fig1_cohort_qc_global_profile"), 183, 135)

# Figure 2: differential landscape.
p_vol1 <- make_volcano(diff_res$E_vs_H, "T1-T2 GC vs healthy")
p_vol2 <- make_volcano(diff_res$L_vs_E, "T3-T4 vs T1-T2 GC")
p_vol3 <- make_volcano(diff_res$L_vs_H, "T3-T4 GC vs healthy")

sig_count <- bind_rows(lapply(names(diff_res), function(nm) {
  diff_res[[nm]] |>
    mutate(Comparison = nm,
           Direction = case_when(adj.P.Val < 0.05 & logFC > 0 ~ "Up",
                                 adj.P.Val < 0.05 & logFC < 0 ~ "Down",
                                 TRUE ~ "NS")) |>
    count(Comparison, Direction)
})) |>
  filter(Direction != "NS") |>
  mutate(Comparison = recode(Comparison, E_vs_H = "T1-T2 vs H", L_vs_E = "T3-T4 vs T1-T2", L_vs_H = "T3-T4 vs H"))
p_sigcount <- ggplot(sig_count, aes(Comparison, n, fill = Direction)) +
  geom_col(width = 0.68, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c(Up = "#D24B40", Down = "#3182BD")) +
  labs(title = "FDR < 0.05 proteins", x = NULL, y = "Count", fill = NULL) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6.8),
    legend.key.width = unit(5.2, "mm"),
    legend.spacing.x = unit(1.2, "mm"),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

top_heat_proteins <- unique(c(head(early_sig$Protein, 8), head(late_sig$Protein, 8), head(hl_sig$Protein, 10)))
top_heat_proteins <- top_heat_proteins[seq_len(min(24, length(top_heat_proteins)))]
p_diff_heat <- make_heatmap(expr, meta, top_heat_proteins, "Top differential proteins")

fig2 <- ((p_vol1 | p_vol2 | p_vol3) / (p_sigcount | p_diff_heat)) +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = "a")
save_pub(fig2, file.path(fig_dir, "Fig2_differential_protein_landscape"), 183, 150)

# Figure 3: ordered progression signature.
trend_plot_proteins <- trend_sig |>
  arrange(trend_FDR) |>
  slice_head(n = min(12, nrow(trend_sig))) |>
  pull(Protein)
trend_line_df <- as.data.frame(expr[, trend_plot_proteins, drop = FALSE]) |>
  mutate(SampleID = rownames(expr), Group = meta$GroupLabel) |>
  pivot_longer(cols = all_of(trend_plot_proteins), names_to = "Protein", values_to = "NPX") |>
  mutate(Protein = factor(Protein, levels = trend_plot_proteins))
p_trend_line <- ggplot(trend_line_df, aes(Group, NPX, group = Protein, colour = Protein)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.45, alpha = 0.9) +
  stat_summary(fun = mean, geom = "point", size = 1.5) +
  labs(title = "Top monotonic T-category proteins", x = NULL, y = "Mean NPX") +
  theme(legend.position = "right")

trend_top_heat <- trend_sig |>
  arrange(trend_FDR) |>
  slice_head(n = min(30, nrow(trend_sig))) |>
  pull(Protein)
p_trend_heat <- make_heatmap(expr, meta, trend_top_heat, "Monotonic T-category protein signature")

p_trend_vol_df <- trend_res |>
  mutate(Signature = case_when(monotonic & trend_FDR < 0.05 & slope > 0 ~ "Increasing",
                               monotonic & trend_FDR < 0.05 & slope < 0 ~ "Decreasing",
                               TRUE ~ "Other"),
         label = ifelse(Signature != "Other" & rank(trend_FDR, ties.method = "first") <= 12, Protein, NA))
p_trend_vol_label <- p_trend_vol_df |> filter(!is.na(label))
p_trend_vol <- p_trend_vol_df |>
  ggplot(aes(slope, -log10(trend_FDR), colour = Signature)) +
  geom_point(size = 0.9, alpha = 0.82) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.25, colour = palette["neutral_mid"]) +
  ggrepel::geom_text_repel(data = p_trend_vol_label, aes(label = label), size = 2, min.segment.length = 0, max.overlaps = Inf) +
  scale_colour_manual(values = c(Increasing = "#D24B40", Decreasing = "#3182BD", Other = "#D8D8D8")) +
  labs(title = "Ordered T-category trend", x = "Slope per group step", y = "-log10 trend FDR") +
  theme(legend.position = "bottom")

fig3 <- ((p_trend_vol | p_trend_line) / p_trend_heat) +
  plot_layout(heights = c(0.9, 1.25)) +
  plot_annotation(tag_levels = "a")
save_pub(fig3, file.path(fig_dir, "Fig3_stage_progression_signature"), 183, 145)

# Figure 4: internally validated sparse diagnostic models.
roc_df <- bind_rows(plot_roc_df(model_early), plot_roc_df(model_late))
roc_labels <- c(
  `T1-T2 GC vs healthy` = "T1-T2 model",
  `T3-T4 vs T1-T2 GC` = "Invasion-depth model"
)
p_roc <- ggplot(roc_df, aes(FPR, Sensitivity, colour = Model)) +
  geom_line(linewidth = 0.65) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.25, colour = palette["neutral_mid"]) +
  coord_equal() +
  scale_colour_manual(values = model_colors, labels = roc_labels) +
  labs(title = "Nested LOOCV performance", x = "1 - specificity", y = "Sensitivity", colour = NULL) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 0.8))) +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.text = element_text(size = 5.8, lineheight = 0.9),
    legend.key.width = unit(6.5, "mm"),
    legend.box.margin = margin(t = 1, r = 1, b = 1, l = 0, unit = "mm")
  )

pred_df <- bind_rows(
  data.frame(SampleID = model_early$meta$SampleID, Group = model_early$meta$GroupLabel, Probability = model_early$pred_cv, Model = model_early$label),
  data.frame(SampleID = model_late$meta$SampleID, Group = model_late$meta$GroupLabel, Probability = model_late$pred_cv, Model = model_late$label)
)
p_pred <- ggplot(pred_df, aes(Group, Probability, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.size = 0.6, colour = "black", linewidth = 0.25) +
  geom_jitter(width = 0.08, size = 0.8, alpha = 0.75) +
  facet_wrap(~ Model, scales = "free_x") +
  scale_fill_manual(values = group_colors) +
  labs(title = "Cross-validated predicted probability", x = NULL, y = "Probability of positive class") +
  theme(legend.position = "none")

freq_df <- bind_rows(
  model_early$feature_frequency |> mutate(Model = model_early$label),
  model_late$feature_frequency |> mutate(Model = model_late$label)
) |>
  group_by(Model) |>
  slice_max(Frequency, n = 12, with_ties = FALSE) |>
  ungroup()
p_freq <- ggplot(freq_df, aes(reorder(Protein, Frequency), Frequency, fill = Model)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.15) +
  coord_flip() +
  facet_wrap(~ Model, scales = "free_y") +
  scale_y_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = model_colors) +
  labs(title = "Feature-selection stability", x = NULL, y = "LOOCV selection frequency") +
  theme(legend.position = "none")

fig4 <- ((p_roc | p_pred) / p_freq) +
  plot_layout(heights = c(1.05, 1.1), widths = c(1.05, 1.75)) +
  plot_annotation(tag_levels = "a")
save_pub(fig4, file.path(fig_dir, "Fig4_sparse_diagnostic_models"), 210, 154)

# Figure 5: pathway context.
if (nrow(go_res) > 0) {
  go_plot <- go_res |>
    group_by(GeneSet) |>
    slice_min(p.adjust, n = 8, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      DescriptionWrapped = stringr::str_wrap(Description, width = 34),
      GeneSetLabel = recode(GeneSet, Early_vs_Healthy = "T1-T2 vs healthy", Late_vs_Early = "T3-T4 vs T1-T2", Late_vs_Healthy = "T3-T4 vs healthy", Monotonic_trend = "Monotonic trend"),
      DescriptionWrapped = reorder(DescriptionWrapped, -log10(p.adjust))
    )
  p_go <- ggplot(go_plot, aes(-log10(p.adjust), DescriptionWrapped, size = Count, colour = GeneSetLabel)) +
    geom_point(alpha = 0.9) +
    facet_wrap(~ GeneSetLabel, scales = "free_y", ncol = 2, labeller = labeller(GeneSetLabel = label_wrap_gen(width = 14))) +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.22))) +
    scale_colour_manual(values = c(`T1-T2 vs healthy` = "#3182BD", `T3-T4 vs T1-T2` = "#D24B40", `T3-T4 vs healthy` = "#E28E2C", `Monotonic trend` = "#33B5A5")) +
    labs(title = "GO biological-process enrichment", x = "-log10 adjusted P", y = NULL, size = "Genes") +
    guides(colour = "none", size = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      legend.box.margin = margin(t = 1, r = 8, b = 1, l = 1, unit = "mm"),
      strip.text = element_text(size = 6.0, lineheight = 0.9, margin = margin(t = 1.5, r = 2, b = 1.5, l = 2, unit = "mm")),
      axis.text.y = element_text(size = 5.6, lineheight = 0.9, margin = margin(r = 3)),
      plot.margin = margin(t = 3, r = 10, b = 4, l = 14, unit = "mm")
    )
  save_pub(p_go, file.path(fig_dir, "Fig5_go_biological_context"), 230, 160)
}

report_lines <- c(
  "# Olink Reveal primary bioinformatics analysis report",
  "",
  "## Cohort and QC",
  paste0("- Samples: ", nrow(expr), " (", paste(names(table(meta$Group)), as.integer(table(meta$Group)), sep = "=", collapse = "; "), ")."),
  paste0("- Protein assays analysed: ", ncol(expr), "."),
  paste0("- Missing NPX values after reading provided matrix: ", sum(is.na(expr_raw)), ". Median imputation was applied only where missing values existed."),
  "",
  "## Differential proteins",
  paste0("- T1-T2 GC vs healthy: ", nrow(early_sig), " proteins at FDR < 0.05."),
  paste0("- T3-T4 vs T1-T2 GC: ", nrow(late_sig), " proteins at FDR < 0.05."),
  paste0("- T3-T4 GC vs healthy: ", nrow(hl_sig), " proteins at FDR < 0.05."),
  "",
  "## Stage-trend proteins",
  paste0("- Monotonic T-category proteins at trend FDR < 0.05: ", nrow(trend_sig), "."),
  paste0("- Top trend proteins: ", paste(head(trend_sig$Protein, 20), collapse = ", "), "."),
  "",
  "## Sparse internally validated models",
  paste0("- ", model_early$label, ": ", auc_ci_text(model_early$roc), ". Final features: ", paste(setdiff(model_early$coef$Protein, "(Intercept)"), collapse = ", "), "."),
  paste0("- ", model_late$label, ": ", auc_ci_text(model_late$roc), ". Final features: ", paste(setdiff(model_late$coef$Protein, "(Intercept)"), collapse = ", "), "."),
  "",
  "## Enrichment",
  paste0("- GO biological-process enrichment rows exported: ", nrow(go_res), "."),
  "",
  "## Review-sensitive interpretation",
  "- These are candidate plasma protein signatures from a 90-sample discovery cohort. Clinical diagnostic utility requires independent validation.",
  "- The nested LOOCV estimates reduce feature-selection leakage, but they do not replace an external validation cohort.",
  "- Mechanistic claims should be phrased as pathway-context hypotheses unless supported by experimental validation."
)
writeLines(report_lines, file.path(report_dir, "primary_analysis_report.md"), useBytes = TRUE)

message("Primary Olink analysis complete.")
message("Tables: ", tab_dir)
message("Figures: ", fig_dir)
message("Report: ", file.path(report_dir, "primary_analysis_report.md"))
