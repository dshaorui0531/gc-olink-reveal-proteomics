#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(pROC)
  library(survival)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(curl)
  library(ragg)
  library(svglite)
})

set.seed(20260819)
options(stringsAsFactors = FALSE)
figure_font <- "Helvetica"

out_dir <- Sys.getenv("GC_OLINK_RESULTS_DIR", "results")
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
report_dir <- file.path(out_dir, "reports")
cache_dir <- file.path(out_dir, "cache", "tcga_stad")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = file.path(out_dir, "cache"))

theme_pub <- function(base_size = 7, base_family = figure_font) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.6),
      strip.text = element_text(size = base_size - 0.2, face = "bold"),
      plot.title = element_text(size = base_size + 0.6, face = "bold"),
      plot.margin = margin(t = 3, r = 7, b = 3, l = 3, unit = "mm"),
      panel.grid = element_blank()
    )
}
theme_set(theme_pub())

save_pub <- function(plot, stem, width_mm = 183, height_mm = 130, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  message("Saving figure: ", basename(stem))
  grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h, family = figure_font, useDingbats = FALSE)
  print(plot)
  grDevices::dev.off()
  ggplot2::ggsave(paste0(stem, ".tiff"), plot, width = w, height = h, dpi = dpi, device = ragg::agg_tiff)
}

download_if_missing <- function(url, dest) {
  if (file.exists(dest) && file.info(dest)$size > 0) {
    message("Using cached file: ", basename(dest))
    return(invisible(dest))
  }
  message("Downloading: ", url)
  curl::curl_download(url, destfile = dest, quiet = FALSE)
  invisible(dest)
}

safe_write <- function(x, file) {
  write.csv(x, file, row.names = FALSE, fileEncoding = "UTF-8")
}

normalise_stage <- function(x) {
  y <- toupper(as.character(x))
  y <- gsub("\\[|\\]|'|\"", "", y)
  y <- trimws(y)
  case_when(
    grepl("STAGE IV", y) ~ "Stage IV",
    grepl("STAGE III", y) ~ "Stage III",
    grepl("STAGE II", y) ~ "Stage II",
    grepl("STAGE I", y) ~ "Stage I",
    TRUE ~ NA_character_
  )
}

expr_url <- "https://gdc.xenahubs.net/download/TCGA-STAD.star_tpm.tsv.gz"
clinical_url <- "https://gdc.xenahubs.net/download/TCGA-STAD.clinical.tsv.gz"
survival_url <- "https://gdc.xenahubs.net/download/TCGA-STAD.survival.tsv.gz"

expr_file <- file.path(cache_dir, "TCGA-STAD.star_tpm.tsv.gz")
clinical_file <- file.path(cache_dir, "TCGA-STAD.clinical.tsv.gz")
survival_file <- file.path(cache_dir, "TCGA-STAD.survival.tsv.gz")

download_if_missing(expr_url, expr_file)
download_if_missing(clinical_url, clinical_file)
download_if_missing(survival_url, survival_file)

early_model <- read.csv(file.path(tab_dir, "model_early_vs_healthy_final_coefficients.csv"))
late_model <- read.csv(file.path(tab_dir, "model_late_vs_early_final_coefficients.csv"))
trend <- read.csv(file.path(tab_dir, "stage_trend_analysis.csv"))

candidate_genes <- unique(c(
  setdiff(early_model$Protein, "(Intercept)"),
  setdiff(late_model$Protein, "(Intercept)"),
  head(trend$Protein[trend$monotonic & trend$trend_FDR < 0.05], 25)
))

message("Reading TCGA-STAD STAR TPM matrix.")
expr_dt <- data.table::fread(expr_file, data.table = FALSE, check.names = FALSE)
ensembl <- sub("\\..*$", "", expr_dt[[1]])
symbols <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = ensembl,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)
keep_rows <- !is.na(symbols) & symbols %in% candidate_genes
expr_sub <- expr_dt[keep_rows, , drop = FALSE]
expr_sub$Symbol <- unname(symbols[keep_rows])
expr_sub[[1]] <- NULL

expr_long <- expr_sub |>
  pivot_longer(cols = -Symbol, names_to = "Sample", values_to = "TPM_log2") |>
  mutate(
    TPM_log2 = as.numeric(TPM_log2),
    SampleTypeCode = substr(Sample, 14, 15),
    Tissue = case_when(
      SampleTypeCode == "01" ~ "Tumor",
      SampleTypeCode == "11" ~ "Adjacent normal",
      TRUE ~ "Other"
    ),
    Patient = substr(Sample, 1, 12)
  ) |>
  filter(Tissue %in% c("Tumor", "Adjacent normal"))

expr_long <- expr_long |>
  group_by(Symbol, Sample, Tissue, Patient) |>
  summarise(TPM_log2 = mean(TPM_log2, na.rm = TRUE), .groups = "drop")

safe_write(expr_long, file.path(tab_dir, "tcga_stad_candidate_expression_long.csv"))

tumor_normal <- expr_long |>
  filter(Tissue %in% c("Tumor", "Adjacent normal")) |>
  group_by(Symbol) |>
  summarise(
    n_tumor = sum(Tissue == "Tumor"),
    n_normal = sum(Tissue == "Adjacent normal"),
    mean_tumor = mean(TPM_log2[Tissue == "Tumor"], na.rm = TRUE),
    mean_normal = mean(TPM_log2[Tissue == "Adjacent normal"], na.rm = TRUE),
    delta_tumor_minus_normal = mean_tumor - mean_normal,
    wilcox_p = ifelse(n_tumor > 2 && n_normal > 2, wilcox.test(TPM_log2 ~ Tissue)$p.value, NA_real_),
    .groups = "drop"
  ) |>
  mutate(FDR = p.adjust(wilcox_p, method = "BH")) |>
  arrange(FDR)
safe_write(tumor_normal, file.path(tab_dir, "tcga_stad_tumor_vs_normal_validation.csv"))

clinical <- data.table::fread(clinical_file, data.table = FALSE, check.names = FALSE)
stage_cols <- grep("ajcc.*stage|pathologic.*stage|tumor_stage", names(clinical), ignore.case = TRUE, value = TRUE)
stage_col <- stage_cols[1]
if (is.na(stage_col)) {
  clinical_stage <- clinical |> transmute(Sample = sample, Patient = substr(sample, 1, 12), Stage = NA_character_)
} else {
  clinical_stage <- clinical |>
    transmute(Sample = sample, Patient = substr(sample, 1, 12), Stage = normalise_stage(.data[[stage_col]])) |>
    filter(!is.na(Stage)) |>
    distinct(Patient, .keep_all = TRUE)
}

stage_expr <- expr_long |>
  filter(Tissue == "Tumor") |>
  inner_join(clinical_stage |> dplyr::select(Patient, Stage), by = "Patient") |>
  mutate(StageOrdinal = as.numeric(factor(Stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))))
safe_write(stage_expr, file.path(tab_dir, "tcga_stad_stage_expression_long.csv"))

stage_assoc <- stage_expr |>
  group_by(Symbol) |>
  summarise(
    n_stage = sum(!is.na(StageOrdinal)),
    spearman_rho = ifelse(n_stage >= 20, suppressWarnings(cor(TPM_log2, StageOrdinal, method = "spearman", use = "complete.obs")), NA_real_),
    spearman_p = ifelse(n_stage >= 20, suppressWarnings(cor.test(TPM_log2, StageOrdinal, method = "spearman")$p.value), NA_real_),
    kruskal_p = ifelse(length(unique(Stage[!is.na(Stage)])) >= 2, kruskal.test(TPM_log2 ~ Stage)$p.value, NA_real_),
    .groups = "drop"
  ) |>
  mutate(spearman_FDR = p.adjust(spearman_p, method = "BH"),
         kruskal_FDR = p.adjust(kruskal_p, method = "BH")) |>
  arrange(spearman_FDR)
safe_write(stage_assoc, file.path(tab_dir, "tcga_stad_stage_association.csv"))

surv <- data.table::fread(survival_file, data.table = FALSE, check.names = FALSE)
surv <- surv |>
  transmute(Sample = sample, Patient = `_PATIENT`, OS_time = as.numeric(OS.time), OS = as.numeric(OS)) |>
  filter(!is.na(OS_time), !is.na(OS), OS_time > 30) |>
  arrange(Patient, desc(OS_time)) |>
  distinct(Patient, .keep_all = TRUE)

tumor_expr_wide <- expr_long |>
  filter(Tissue == "Tumor") |>
  group_by(Patient, Symbol) |>
  summarise(TPM_log2 = mean(TPM_log2, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = Symbol, values_from = TPM_log2)

surv_expr <- surv |>
  inner_join(tumor_expr_wide, by = "Patient")

surv_gene_rows <- lapply(intersect(candidate_genes, names(surv_expr)), function(g) {
  df <- surv_expr |> filter(!is.na(.data[[g]]))
  if (nrow(df) < 40 || length(unique(df$OS)) < 2) return(NULL)
  med <- median(df[[g]], na.rm = TRUE)
  df$High <- as.integer(df[[g]] > med)
  cox <- survival::coxph(Surv(OS_time, OS) ~ High, data = df)
  s <- summary(cox)
  data.frame(
    Gene = g,
    n = nrow(df),
    HR_high_vs_low = unname(exp(coef(cox))),
    HR_lower95 = s$conf.int[1, "lower .95"],
    HR_upper95 = s$conf.int[1, "upper .95"],
    cox_p = s$coefficients[1, "Pr(>|z|)"],
    row.names = NULL
  )
}) |> bind_rows()
if (nrow(surv_gene_rows) > 0) {
  surv_gene_rows$FDR <- p.adjust(surv_gene_rows$cox_p, method = "BH")
  surv_gene_rows <- surv_gene_rows |> arrange(cox_p)
}
safe_write(surv_gene_rows, file.path(tab_dir, "tcga_stad_candidate_survival.csv"))

make_signature_score <- function(df, coef_table, score_name) {
  coef_table <- coef_table |> filter(Protein != "(Intercept)")
  genes <- intersect(coef_table$Protein, names(df))
  if (length(genes) < 2) return(NULL)
  x <- as.matrix(df[, genes, drop = FALSE])
  x <- scale(x)
  beta <- coef_table$Coefficient[match(genes, coef_table$Protein)]
  score <- as.numeric(x %*% beta)
  data.frame(Patient = df$Patient, Score = score, Signature = score_name)
}

sig_scores <- bind_rows(
  make_signature_score(surv_expr, early_model, "T1-T2 detection signature"),
  make_signature_score(surv_expr, late_model, "Invasion-depth signature")
) |>
  inner_join(surv |> dplyr::select(Patient, OS_time, OS), by = "Patient")

sig_surv <- sig_scores |>
  group_by(Signature) |>
  group_modify(~ {
    med <- median(.x$Score, na.rm = TRUE)
    .x$High <- as.integer(.x$Score > med)
    cox <- survival::coxph(Surv(OS_time, OS) ~ High, data = .x)
    s <- summary(cox)
    data.frame(
      n = nrow(.x),
      HR_high_vs_low = unname(exp(coef(cox))),
      HR_lower95 = s$conf.int[1, "lower .95"],
      HR_upper95 = s$conf.int[1, "upper .95"],
      cox_p = s$coefficients[1, "Pr(>|z|)"]
    )
  }) |>
  ungroup() |>
  mutate(FDR = p.adjust(cox_p, method = "BH"))
safe_write(sig_surv, file.path(tab_dir, "tcga_stad_signature_survival.csv"))
safe_write(sig_scores, file.path(tab_dir, "tcga_stad_signature_scores.csv"))

plot_genes <- unique(c(setdiff(early_model$Protein, "(Intercept)"), setdiff(late_model$Protein, "(Intercept)")))
p_tn <- expr_long |>
  filter(Symbol %in% plot_genes, Tissue %in% c("Adjacent normal", "Tumor")) |>
  ggplot(aes(Tissue, TPM_log2, fill = Tissue)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3, colour = "black", linewidth = 0.25) +
  facet_wrap(~ Symbol, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c("Adjacent normal" = "#D8D8D8", Tumor = "#D24B40")) +
  labs(title = "TCGA-STAD tumor versus adjacent normal tissue", x = NULL, y = "STAR TPM (log2)") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 25, hjust = 1), plot.margin = margin(t = 3, r = 9, b = 3, l = 3, unit = "mm"))

stage_plot_candidates <- stage_assoc |>
  filter(!is.na(spearman_p)) |>
  arrange(spearman_p)
stage_plot_genes <- stage_plot_candidates |>
  slice_head(n = min(8, nrow(stage_plot_candidates))) |>
  pull(Symbol)
p_stage <- stage_expr |>
  filter(Symbol %in% stage_plot_genes) |>
  ggplot(aes(Stage, TPM_log2, fill = Stage)) +
  geom_boxplot(width = 0.6, outlier.size = 0.3, colour = "black", linewidth = 0.25) +
  facet_wrap(~ Symbol, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Stage I" = "#D8D8D8", "Stage II" = "#9ECAE1", "Stage III" = "#3182BD", "Stage IV" = "#D24B40")) +
  labs(title = "Candidate gene expression across pathological stage", x = NULL, y = "STAR TPM (log2)") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1), plot.margin = margin(t = 3, r = 9, b = 3, l = 3, unit = "mm"))

forest_candidates <- surv_gene_rows |>
  filter(!is.na(cox_p)) |>
  arrange(cox_p)
forest_df <- forest_candidates |>
  slice_head(n = min(14, nrow(forest_candidates))) |>
  mutate(Gene = factor(Gene, levels = rev(Gene)))
p_forest <- ggplot(forest_df, aes(HR_high_vs_low, Gene)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#767676", linewidth = 0.25) +
  geom_errorbar(aes(xmin = HR_lower95, xmax = HR_upper95), orientation = "y", width = 0.18, linewidth = 0.35) +
  geom_point(size = 1.6, colour = "#272727") +
  scale_x_log10(expand = expansion(mult = c(0.04, 0.18))) +
  labs(title = "Overall survival association", x = "Hazard ratio, high vs low expression", y = NULL)

p_sig <- sig_scores |>
  group_by(Signature) |>
  mutate(Group = ifelse(Score > median(Score, na.rm = TRUE), "High", "Low")) |>
  ungroup() |>
  ggplot(aes(Group, Score, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.size = 0.5, colour = "black", linewidth = 0.25) +
  facet_wrap(~ Signature, scales = "free_y", labeller = labeller(Signature = label_wrap_gen(width = 24))) +
  scale_fill_manual(values = c(High = "#D24B40", Low = "#3182BD")) +
  labs(title = "TCGA mRNA signature scores", x = NULL, y = "Coefficient-weighted z score") +
  theme(legend.position = "none", plot.margin = margin(t = 3, r = 10, b = 3, l = 3, unit = "mm"))

fig6 <- (p_tn / p_stage) / (p_forest | p_sig) +
  plot_layout(heights = c(1.1, 1.0, 0.9)) +
  plot_annotation(tag_levels = "a")
save_pub(fig6, file.path(fig_dir, "Fig6_tcga_stad_external_validation"), 205, 218)

report_lines <- c(
  "# TCGA-STAD external validation report",
  "",
  "## Data source",
  "- UCSC Xena/GDC TCGA-STAD STAR TPM expression matrix, clinical table, and survival table were downloaded from the current GDC Xena hub.",
  paste0("- Candidate genes tested: ", paste(candidate_genes, collapse = ", "), "."),
  "",
  "## Tumor versus adjacent normal",
  paste0("- Candidate genes with tumor-normal FDR < 0.05: ", sum(tumor_normal$FDR < 0.05, na.rm = TRUE), " of ", nrow(tumor_normal), "."),
  paste0("- Top tumor-normal genes: ", paste(head(tumor_normal$Symbol, 10), collapse = ", "), "."),
  "",
  "## Pathological stage",
  paste0("- Stage column used: ", ifelse(is.na(stage_col), "none found", stage_col), "."),
  paste0("- Candidate genes with Spearman stage FDR < 0.05: ", sum(stage_assoc$spearman_FDR < 0.05, na.rm = TRUE), " of ", nrow(stage_assoc), "."),
  "",
  "## Overall survival",
  paste0("- Candidate genes with Cox FDR < 0.05: ", sum(surv_gene_rows$FDR < 0.05, na.rm = TRUE), " of ", nrow(surv_gene_rows), "."),
  paste0("- Signature survival rows: ", nrow(sig_surv), "."),
  "",
  "## Interpretation boundary",
  "- This validation is tissue mRNA-level support for candidate genes; it does not validate plasma protein concentrations.",
  "- Survival and stage analyses are exploratory and should be presented as orthogonal public-dataset context."
)
writeLines(report_lines, file.path(report_dir, "tcga_stad_external_validation_report.md"), useBytes = TRUE)

message("TCGA-STAD external validation complete.")
message("Report: ", file.path(report_dir, "tcga_stad_external_validation_report.md"))
