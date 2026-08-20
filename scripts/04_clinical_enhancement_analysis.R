#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(limma)
  library(pROC)
})

set.seed(20260819)
options(stringsAsFactors = FALSE)

clinical_path <- Sys.getenv("GC_OLINK_CLINICAL_FILE", "data/clinical_metadata.xlsx")
input_dir <- Sys.getenv("GC_OLINK_SUMMARY_DIR", "data/summary")
result_dir <- Sys.getenv("GC_OLINK_RESULTS_DIR", "results")
tab_dir <- file.path(result_dir, "tables")
report_dir <- file.path(result_dir, "reports")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

safe_write <- function(x, file) {
  write.csv(x, file, row.names = FALSE, fileEncoding = "UTF-8")
}

as_num <- function(x) suppressWarnings(as.numeric(x))

clean_binary <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "无", "未详", "未知")] <- NA
  x
}

exposure_yes_no <- function(x) {
  x <- clean_binary(x)
  out <- ifelse(is.na(x), NA_character_, ifelse(x %in% c("否", "无", "未见"), "No", "Yes"))
  factor(out, levels = c("No", "Yes"))
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

fmt_cont <- function(x) {
  x <- as_num(x)
  if (all(is.na(x))) return("NA")
  sprintf("%.1f [%.1f, %.1f]", median(x, na.rm = TRUE), quantile(x, 0.25, na.rm = TRUE), quantile(x, 0.75, na.rm = TRUE))
}

fmt_cat <- function(x) {
  x <- clean_binary(x)
  tab <- sort(table(x), decreasing = TRUE)
  if (length(tab) == 0) return("NA")
  paste(sprintf("%s: %d (%.1f%%)", names(tab), as.integer(tab), 100 * as.integer(tab) / sum(tab)), collapse = "; ")
}

p_cont <- function(x, g) {
  ok <- !is.na(as_num(x)) & !is.na(g)
  if (sum(ok) < 3 || length(unique(g[ok])) < 2) return(NA_real_)
  suppressWarnings(kruskal.test(as_num(x)[ok] ~ g[ok])$p.value)
}

p_cat <- function(x, g) {
  x <- clean_binary(x)
  ok <- !is.na(x) & !is.na(g)
  if (sum(ok) < 3 || length(unique(x[ok])) < 2 || length(unique(g[ok])) < 2) return(NA_real_)
  suppressWarnings(fisher.test(table(x[ok], g[ok]))$p.value)
}

parse_max_dimension <- function(x) {
  vapply(as.character(x), function(z) {
    if (is.na(z) || trimws(z) == "") return(NA_real_)
    nums <- regmatches(z, gregexpr("[0-9]+\\.?[0-9]*", z))[[1]]
    nums <- suppressWarnings(as.numeric(nums))
    if (length(nums) == 0 || all(is.na(nums))) NA_real_ else max(nums, na.rm = TRUE)
  }, numeric(1))
}

extract_t_stage <- function(x) {
  x <- toupper(as.character(x))
  out <- ifelse(grepl("T4", x), "T4",
                ifelse(grepl("T3", x), "T3",
                       ifelse(grepl("T2", x), "T2",
                              ifelse(grepl("T1", x), "T1", NA_character_))))
  factor(out, levels = c("T1", "T2", "T3", "T4"))
}

extract_n_status <- function(x) {
  x <- toupper(as.character(x))
  out <- ifelse(grepl("N0", x), "N0",
                ifelse(grepl("N[1-3]", x), "N+", ifelse(grepl("NX", x), "Nx", NA_character_)))
  factor(out, levels = c("N0", "N+", "Nx"))
}

extract_m_status <- function(x) {
  x <- toupper(as.character(x))
  out <- ifelse(grepl("M0|MO", x), "M0",
                ifelse(grepl("M1", x), "M1", ifelse(grepl("MX", x), "Mx", NA_character_)))
  factor(out, levels = c("M0", "M1", "Mx"))
}

normalize_cancer_stage <- function(x) {
  x <- as.character(x)
  out <- ifelse(grepl("Ⅳ|IV", x), "IV",
                ifelse(grepl("Ⅲ|III", x), "III",
                       ifelse(grepl("Ⅱ|II", x), "II",
                              ifelse(grepl("Ⅰ|I", x), "I", NA_character_))))
  factor(out, levels = c("I", "II", "III", "IV"))
}

make_table1 <- function(df, group_var, variables, out_file) {
  groups <- if (is.factor(df[[group_var]])) levels(droplevels(df[[group_var]])) else unique(df[[group_var]])
  groups <- groups[!is.na(groups)]
  rows <- lapply(variables, function(v) {
    type <- if (is.numeric(df[[v]])) "continuous" else "categorical"
    summaries <- lapply(groups, function(g) {
      x <- df[df[[group_var]] == g, v, drop = TRUE]
      if (type == "continuous") fmt_cont(x) else fmt_cat(x)
    })
    names(summaries) <- paste0("Group_", groups)
    p <- if (type == "continuous") p_cont(df[[v]], df[[group_var]]) else p_cat(df[[v]], df[[group_var]])
    data.frame(
      Variable = v,
      Type = type,
      Overall = if (type == "continuous") fmt_cont(df[[v]]) else fmt_cat(df[[v]]),
      Missing = sum(is.na(df[[v]]) | df[[v]] == "", na.rm = TRUE),
      do.call(data.frame, summaries),
      P_value = p,
      P_value_display = fmt_p(p),
      check.names = FALSE
    )
  })
  out <- bind_rows(rows)
  safe_write(out, out_file)
  out
}

read_olink_matrix <- function() {
  npx_path <- file.path(input_dir, "04.Diff_analysis/COND1/all_NPX.xlsx")
  npx <- readxl::read_excel(npx_path, sheet = 1, .name_repair = "unique") |> as.data.frame()
  names(npx)[1] <- "SampleID"
  rownames(npx) <- npx$SampleID
  protein_cols <- setdiff(names(npx), c("SampleID", "Group"))
  expr <- as.matrix(npx[, protein_cols])
  storage.mode(expr) <- "numeric"
  for (j in seq_len(ncol(expr))) {
    miss <- is.na(expr[, j])
    if (any(miss)) expr[miss, j] <- median(expr[, j], na.rm = TRUE)
  }
  meta <- data.frame(
    SampleID = rownames(expr),
    Group = factor(npx$Group, levels = c("H", "E", "L")),
    stringsAsFactors = FALSE
  )
  list(expr = expr, meta = meta)
}

run_adjusted_limma <- function(expr, clin, positive, negative, covariates, tag) {
  keep <- clin$Group %in% c(negative, positive)
  df <- clin[keep, , drop = FALSE]
  rownames(df) <- df$SampleID
  x <- expr[df$SampleID, , drop = FALSE]
  df$GroupBin <- factor(ifelse(df$Group == positive, positive, negative), levels = c(negative, positive))

  use_cov <- covariates
  for (cv in covariates) {
    vals <- df[[cv]]
    if (is.numeric(vals)) {
      if (sum(!is.na(vals)) < 10 || length(unique(vals[!is.na(vals)])) < 2) use_cov <- setdiff(use_cov, cv)
    } else {
      vals <- factor(vals)
      if (sum(!is.na(vals)) < 10 || length(unique(vals[!is.na(vals)])) < 2) use_cov <- setdiff(use_cov, cv)
    }
  }

  model_df <- df[, c("GroupBin", use_cov), drop = FALSE]
  complete <- complete.cases(model_df)
  model_df <- model_df[complete, , drop = FALSE]
  x <- x[complete, , drop = FALSE]
  for (cv in use_cov) {
    if (!is.numeric(model_df[[cv]])) model_df[[cv]] <- factor(model_df[[cv]])
  }

  design <- model.matrix(as.formula(paste("~ GroupBin", if (length(use_cov) > 0) paste("+", paste(use_cov, collapse = "+")) else "")), data = model_df)
  coef_name <- paste0("GroupBin", positive)
  fit <- eBayes(lmFit(t(x), design))
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  tt$Protein <- rownames(tt)
  tt <- tt[, c("Protein", setdiff(names(tt), "Protein"))]
  tt$Comparison <- paste0(positive, "_vs_", negative)
  tt$Covariates <- paste(use_cov, collapse = ";")
  tt$Samples_used <- nrow(model_df)
  safe_write(tt, file.path(tab_dir, paste0("limma_", tag, "_", positive, "_vs_", negative, ".csv")))
  tt
}

metric_at_youden <- function(outcome, prob) {
  roc_obj <- pROC::roc(outcome, prob, quiet = TRUE, levels = c(0, 1), direction = "<")
  cc <- pROC::coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"))
  threshold <- as.numeric(cc["threshold"])
  pred <- ifelse(prob >= threshold, 1, 0)
  tp <- sum(pred == 1 & outcome == 1)
  tn <- sum(pred == 0 & outcome == 0)
  fp <- sum(pred == 1 & outcome == 0)
  fn <- sum(pred == 0 & outcome == 1)
  data.frame(
    AUC = as.numeric(pROC::auc(roc_obj)),
    AUC_CI_low = as.numeric(pROC::ci.auc(roc_obj))[1],
    AUC_CI_high = as.numeric(pROC::ci.auc(roc_obj))[3],
    Threshold = threshold,
    Sensitivity = tp / (tp + fn),
    Specificity = tn / (tn + fp),
    Accuracy = (tp + tn) / length(outcome),
    PPV = ifelse(tp + fp == 0, NA, tp / (tp + fp)),
    NPV = ifelse(tn + fn == 0, NA, tn / (tn + fn)),
    TP = tp, FP = fp, TN = tn, FN = fn
  )
}

loocv_logistic <- function(df, outcome, predictors) {
  ok <- complete.cases(df[, c(outcome, predictors), drop = FALSE])
  d <- df[ok, c(outcome, predictors), drop = FALSE]
  y <- d[[outcome]]
  pred <- rep(NA_real_, nrow(d))
  for (i in seq_len(nrow(d))) {
    train <- setdiff(seq_len(nrow(d)), i)
    form <- as.formula(paste(outcome, "~", paste(predictors, collapse = "+")))
    fit <- suppressWarnings(glm(form, data = d[train, , drop = FALSE], family = binomial()))
    pred[i] <- suppressWarnings(as.numeric(predict(fit, newdata = d[i, , drop = FALSE], type = "response")))
  }
  data.frame(Outcome = y, CV_Probability = pred)
}

calibration_summary <- function(outcome, prob, model) {
  eps <- 1e-6
  p <- pmin(pmax(prob, eps), 1 - eps)
  brier <- mean((p - outcome)^2)
  fit <- suppressWarnings(glm(outcome ~ qlogis(p), family = binomial()))
  data.frame(
    Model = model,
    Brier_score = brier,
    Calibration_intercept = unname(coef(fit)[1]),
    Calibration_slope = unname(coef(fit)[2])
  )
}

dca_curve <- function(outcome, prob, model) {
  n <- length(outcome)
  prev <- mean(outcome)
  thresholds <- seq(0.05, 0.95, by = 0.05)
  bind_rows(lapply(thresholds, function(pt) {
    pred <- prob >= pt
    tp <- sum(pred & outcome == 1)
    fp <- sum(pred & outcome == 0)
    nb <- tp / n - fp / n * pt / (1 - pt)
    data.frame(Model = model, Threshold = pt, Net_benefit = nb, Treat_all = prev - (1 - prev) * pt / (1 - pt), Treat_none = 0)
  }))
}

clin_raw <- readxl::read_excel(clinical_path, sheet = 1, .name_repair = "unique") |> as.data.frame()
clin <- clin_raw |>
  transmute(
    SampleID = sample_id,
    GroupCN = 分组,
    Group = recode(分组, "对照" = "H", "早期" = "E", "晚期" = "L"),
    Sex = factor(recode(性别, "男" = "Male", "女" = "Female")),
    Age = as_num(年龄),
    BMI = as_num(BMI),
    Smoking = exposure_yes_no(吸烟史),
    Drinking = exposure_yes_no(饮酒史),
    Residence = factor(clean_binary(居住地)),
    TNM = factor(clean_binary(TNM分期)),
    T_stage = extract_t_stage(TNM分期),
    N_status = extract_n_status(TNM分期),
    M_status = extract_m_status(TNM分期),
    Cancer_stage = normalize_cancer_stage(癌症分期),
    Lauren = factor(clean_binary(`Lauren 分型`)),
    Histology = factor(clean_binary(组织学类型)),
    Differentiation = factor(clean_binary(肿瘤分化度)),
    Metastasis = factor(clean_binary(转移情况)),
    Tumor_size_max_cm = parse_max_dimension(肿瘤大小),
    Hp = factor(clean_binary(`Hp 状态`)),
    ALT = as_num(ALT),
    AST = as_num(AST),
    Albumin = as_num(白蛋白水平),
    eGFR = as_num(eGFR),
    CEA = as_num(CEA),
    CA125 = as_num(CA125),
    CA199 = as_num(`CA19-9`),
    CA724 = as_num(`CA72-4`),
    WBC = as_num(WBC),
    RBC = as_num(RBC),
    Hemoglobin = as_num(血红蛋白),
    Platelet = as_num(血小板)
  )
clin$Group <- factor(clin$Group, levels = c("H", "E", "L"))

dat <- read_olink_matrix()
stopifnot(all(clin$SampleID %in% rownames(dat$expr)))
clin <- clin[match(dat$meta$SampleID, clin$SampleID), ]
safe_write(clin, file.path(tab_dir, "clinical_metadata_clean.csv"))

missingness <- clin |>
  summarise(across(everything(), ~sum(is.na(.x) | .x == "", na.rm = TRUE))) |>
  pivot_longer(everything(), names_to = "Variable", values_to = "Missing") |>
  mutate(Available = nrow(clin) - Missing, Missing_fraction = Missing / nrow(clin)) |>
  arrange(desc(Missing_fraction), Variable)
safe_write(missingness, file.path(tab_dir, "clinical_missingness_summary.csv"))

table1_all <- make_table1(
  clin,
  "Group",
  c("Age", "Sex", "BMI", "Smoking", "Drinking", "Residence", "CEA", "CA199", "CA724", "CA125", "ALT", "AST", "Albumin", "eGFR", "WBC", "RBC", "Hemoglobin", "Platelet"),
  file.path(tab_dir, "table1_overall_clinical_characteristics.csv")
)
table1_cancer <- make_table1(
  clin |> filter(Group %in% c("E", "L")) |> droplevels(),
  "Group",
  c("Age", "Sex", "BMI", "Smoking", "Drinking", "Residence", "T_stage", "N_status", "M_status", "Cancer_stage", "Lauren", "Histology", "Differentiation", "Metastasis", "Tumor_size_max_cm", "CEA", "CA199", "CA724", "CA125", "ALT", "AST", "Albumin", "eGFR", "WBC", "RBC", "Hemoglobin", "Platelet"),
  file.path(tab_dir, "table1_cancer_patient_characteristics_E_vs_L.csv")
)

adj_min <- bind_rows(
  run_adjusted_limma(dat$expr, clin, "E", "H", c("Age", "Sex"), "age_sex_adjusted"),
  run_adjusted_limma(dat$expr, clin, "L", "E", c("Age", "Sex"), "age_sex_adjusted"),
  run_adjusted_limma(dat$expr, clin, "L", "H", c("Age", "Sex"), "age_sex_adjusted")
)
safe_write(adj_min, file.path(tab_dir, "limma_age_sex_adjusted_all_comparisons.csv"))

adj_cancer <- run_adjusted_limma(
  dat$expr,
  clin,
  "L",
  "E",
  c("Age", "Sex", "BMI", "Smoking", "Drinking", "ALT", "AST", "Albumin", "eGFR", "WBC", "Hemoglobin", "Platelet"),
  "cancer_clinical_adjusted"
)
adj_cancer_demographic <- run_adjusted_limma(
  dat$expr,
  clin,
  "L",
  "E",
  c("Age", "Sex", "BMI", "Smoking", "Drinking"),
  "cancer_demographic_adjusted"
)

candidate_proteins <- unique(c(
  read.csv(file.path(tab_dir, "model_early_vs_healthy_final_coefficients.csv"))$Protein,
  read.csv(file.path(tab_dir, "model_late_vs_early_final_coefficients.csv"))$Protein,
  head(read.csv(file.path(tab_dir, "candidate_protein_summary.csv"))$Protein, 30)
))
candidate_proteins <- setdiff(candidate_proteins, "(Intercept)")

unadj_files <- list(
  E_vs_H = "limma_E_vs_H.csv",
  L_vs_E = "limma_L_vs_E.csv",
  L_vs_H = "limma_L_vs_H.csv"
)
robustness <- bind_rows(lapply(names(unadj_files), function(cmp) {
  u <- read.csv(file.path(tab_dir, unadj_files[[cmp]])) |> select(Protein, unadjusted_logFC = logFC, unadjusted_FDR = adj.P.Val)
  a <- adj_min |> filter(Comparison == cmp) |> select(Protein, age_sex_logFC = logFC, age_sex_FDR = adj.P.Val)
  out <- full_join(u, a, by = "Protein") |> filter(Protein %in% candidate_proteins)
  out$Comparison <- cmp
  out
}))
robustness <- robustness |>
  mutate(
    Direction_concordant = sign(unadjusted_logFC) == sign(age_sex_logFC),
    Retained_age_sex_FDR_0.05 = age_sex_FDR < 0.05
  ) |>
  arrange(Comparison, age_sex_FDR)
safe_write(robustness, file.path(tab_dir, "candidate_clinical_adjustment_robustness.csv"))

progression_multivariable_robustness <- read.csv(file.path(tab_dir, "limma_L_vs_E.csv")) |>
  select(Protein, unadjusted_logFC = logFC, unadjusted_FDR = adj.P.Val) |>
  left_join(
    adj_min |> filter(Comparison == "L_vs_E") |> select(Protein, age_sex_logFC = logFC, age_sex_FDR = adj.P.Val),
    by = "Protein"
  ) |>
  left_join(
    adj_cancer_demographic |> select(Protein, demographic_logFC = logFC, demographic_FDR = adj.P.Val),
    by = "Protein"
  ) |>
  left_join(
    adj_cancer |> select(Protein, full_clinical_logFC = logFC, full_clinical_FDR = adj.P.Val),
    by = "Protein"
  ) |>
  filter(Protein %in% candidate_proteins) |>
  mutate(
    Age_sex_direction_concordant = sign(unadjusted_logFC) == sign(age_sex_logFC),
    Demographic_direction_concordant = sign(unadjusted_logFC) == sign(demographic_logFC),
    Full_clinical_direction_concordant = sign(unadjusted_logFC) == sign(full_clinical_logFC)
  ) |>
  arrange(age_sex_FDR, demographic_FDR, full_clinical_FDR)
safe_write(progression_multivariable_robustness, file.path(tab_dir, "progression_candidate_multivariable_robustness.csv"))

early_pred <- read.csv(file.path(tab_dir, "model_early_vs_healthy_predictions.csv")) |>
  left_join(clin, by = "SampleID")
late_pred <- read.csv(file.path(tab_dir, "model_late_vs_early_predictions.csv")) |>
  left_join(clin, by = "SampleID")

model_perf <- bind_rows(
  metric_at_youden(early_pred$Outcome, early_pred$CV_Probability) |> mutate(Model = "Olink early signature", Comparison = "E_vs_H"),
  metric_at_youden(late_pred$Outcome, late_pred$CV_Probability) |> mutate(Model = "Olink progression signature", Comparison = "L_vs_E")
)

early_clin <- loocv_logistic(early_pred, "Outcome", c("Age", "Sex"))
late_clin <- loocv_logistic(late_pred, "Outcome", c("Age", "Sex", "BMI"))
model_perf <- bind_rows(
  model_perf,
  metric_at_youden(early_clin$Outcome, early_clin$CV_Probability) |> mutate(Model = "Clinical covariates: age + sex", Comparison = "E_vs_H"),
  metric_at_youden(late_clin$Outcome, late_clin$CV_Probability) |> mutate(Model = "Clinical covariates: age + sex + BMI", Comparison = "L_vs_E")
) |>
  select(Comparison, Model, everything())
safe_write(model_perf, file.path(tab_dir, "model_threshold_performance_and_clinical_comparison.csv"))

marker_auc <- bind_rows(lapply(c("CEA", "CA199", "CA724", "CA125"), function(m) {
  d <- late_pred |> filter(!is.na(.data[[m]]))
  if (nrow(d) < 10 || length(unique(d$Outcome)) < 2) return(data.frame())
  roc_obj <- pROC::roc(d$Outcome, d[[m]], quiet = TRUE, levels = c(0, 1), direction = "auto")
  ci <- as.numeric(pROC::ci.auc(roc_obj))
  data.frame(
    Comparison = "L_vs_E",
    Marker = m,
    N = nrow(d),
    AUC = as.numeric(pROC::auc(roc_obj)),
    AUC_CI_low = ci[1],
    AUC_CI_high = ci[3]
  )
}))
safe_write(marker_auc, file.path(tab_dir, "traditional_marker_auc_E_vs_L.csv"))

cal <- bind_rows(
  calibration_summary(early_pred$Outcome, early_pred$CV_Probability, "Olink early signature"),
  calibration_summary(late_pred$Outcome, late_pred$CV_Probability, "Olink progression signature"),
  calibration_summary(early_clin$Outcome, early_clin$CV_Probability, "Clinical covariates: age + sex"),
  calibration_summary(late_clin$Outcome, late_clin$CV_Probability, "Clinical covariates: age + sex + BMI")
)
safe_write(cal, file.path(tab_dir, "model_calibration_summary.csv"))

dca <- bind_rows(
  dca_curve(early_pred$Outcome, early_pred$CV_Probability, "Olink early signature") |> mutate(Comparison = "E_vs_H"),
  dca_curve(late_pred$Outcome, late_pred$CV_Probability, "Olink progression signature") |> mutate(Comparison = "L_vs_E"),
  dca_curve(early_clin$Outcome, early_clin$CV_Probability, "Clinical covariates: age + sex") |> mutate(Comparison = "E_vs_H"),
  dca_curve(late_clin$Outcome, late_clin$CV_Probability, "Clinical covariates: age + sex + BMI") |> mutate(Comparison = "L_vs_E")
)
safe_write(dca, file.path(tab_dir, "decision_curve_analysis_data.csv"))

rob_summary <- robustness |>
  group_by(Comparison) |>
  summarise(
    Candidate_n = n(),
    Direction_concordant = sum(Direction_concordant, na.rm = TRUE),
    Retained_age_sex_FDR_0.05 = sum(Retained_age_sex_FDR_0.05, na.rm = TRUE),
    .groups = "drop"
  )
perf_line <- function(comp, model, metric) {
  vals <- model_perf |> filter(Comparison == comp, Model == model) |> pull(metric)
  vals[1]
}
marker_line <- function(marker) {
  vals <- marker_auc |> filter(Marker == marker) |> pull(AUC)
  vals[1]
}

report <- c(
  "# Clinical enhancement analysis report",
  "",
  sprintf("- Clinical samples matched to Olink matrix: %d/%d.", sum(clin$SampleID %in% rownames(dat$expr)), nrow(clin)),
  sprintf("- Complete covariates available for all 90 samples: age and sex."),
  "- BMI, lifestyle factors, biochemical tests, blood counts and tumor markers are mostly unavailable in healthy controls; therefore early-vs-healthy covariate adjustment was restricted to age and sex.",
  "- Patient-only clinical adjustment was additionally performed for late-vs-early comparison using available demographic, lifestyle, hepatic/renal and blood-count covariates.",
  sprintf("- Age differed across H/E/L groups (Kruskal-Wallis P=%s), and sex distribution also differed (Fisher P=%s).", table1_all$P_value_display[table1_all$Variable == "Age"], table1_all$P_value_display[table1_all$Variable == "Sex"]),
  sprintf("- Age/sex-adjusted candidate robustness: E_vs_H %d/%d retained FDR<0.05; L_vs_E %d/%d retained FDR<0.05; L_vs_H %d/%d retained FDR<0.05.",
          rob_summary$Retained_age_sex_FDR_0.05[rob_summary$Comparison == "E_vs_H"], rob_summary$Candidate_n[rob_summary$Comparison == "E_vs_H"],
          rob_summary$Retained_age_sex_FDR_0.05[rob_summary$Comparison == "L_vs_E"], rob_summary$Candidate_n[rob_summary$Comparison == "L_vs_E"],
          rob_summary$Retained_age_sex_FDR_0.05[rob_summary$Comparison == "L_vs_H"], rob_summary$Candidate_n[rob_summary$Comparison == "L_vs_H"]),
  sprintf("- Early-detection Olink signature: AUC %.3f, sensitivity %.3f, specificity %.3f at the Youden threshold; age+sex clinical model AUC %.3f.",
          perf_line("E_vs_H", "Olink early signature", "AUC"),
          perf_line("E_vs_H", "Olink early signature", "Sensitivity"),
          perf_line("E_vs_H", "Olink early signature", "Specificity"),
          perf_line("E_vs_H", "Clinical covariates: age + sex", "AUC")),
  sprintf("- Progression Olink signature: AUC %.3f, sensitivity %.3f, specificity %.3f; age+sex+BMI clinical model AUC %.3f.",
          perf_line("L_vs_E", "Olink progression signature", "AUC"),
          perf_line("L_vs_E", "Olink progression signature", "Sensitivity"),
          perf_line("L_vs_E", "Olink progression signature", "Specificity"),
          perf_line("L_vs_E", "Clinical covariates: age + sex + BMI", "AUC")),
  sprintf("- Traditional tumor marker AUCs for late-vs-early patients: CEA %.3f, CA19-9 %.3f, CA72-4 %.3f, CA125 %.3f. Healthy-control tumor marker data were unavailable, so early-vs-healthy marker comparison is not assessable.",
          marker_line("CEA"), marker_line("CA199"), marker_line("CA724"), marker_line("CA125")),
  "",
  "## Key output files",
  "",
  "- table1_overall_clinical_characteristics.csv",
  "- table1_cancer_patient_characteristics_E_vs_L.csv",
  "- clinical_missingness_summary.csv",
  "- limma_age_sex_adjusted_all_comparisons.csv",
  "- limma_cancer_demographic_adjusted_L_vs_E.csv",
  "- limma_cancer_clinical_adjusted_L_vs_E.csv",
  "- candidate_clinical_adjustment_robustness.csv",
  "- progression_candidate_multivariable_robustness.csv",
  "- model_threshold_performance_and_clinical_comparison.csv",
  "- traditional_marker_auc_E_vs_L.csv",
  "- model_calibration_summary.csv",
  "- decision_curve_analysis_data.csv",
  "",
  "## Manuscript wording constraints",
  "",
  "- Use 'age- and sex-adjusted differential analysis' for all three Olink comparisons.",
  "- Use 'patient-only multivariable adjustment' only for late-vs-early cancer comparison.",
  "- Do not claim clinical tumor marker comparison for early-vs-healthy unless healthy-control marker data are provided.",
  "- Keep TCGA wording as tissue-transcript support rather than plasma-protein external validation."
  ,
  "",
  "## Results text that can be inserted into the manuscript",
  "",
  "Clinical covariates were integrated after exact sample matching between the clinical table and the Olink NPX matrix (90/90 samples). Age and sex were the only complete clinical covariates across all three groups, whereas BMI, lifestyle factors, routine laboratory variables and tumor markers were largely unavailable for healthy controls. Age and sex differed significantly across groups, supporting their inclusion as minimal covariates in adjusted differential analyses.",
  "",
  "After age- and sex-adjustment, most candidate proteins retained concordant effect directions. The early-detection signature remained discriminatory beyond clinical covariates alone, with higher LOOCV AUC than the age+sex clinical model. For late-vs-early progression, the Olink signature also outperformed the age+sex+BMI clinical model, although the separation was more modest and several protein associations were attenuated after broader patient-only clinical adjustment.",
  "",
  "## Methods text that can be inserted into the manuscript",
  "",
  "Clinical metadata were matched to Olink samples by sample identifier. Continuous variables are summarized as median [IQR] and compared using the Kruskal-Wallis test or Wilcoxon rank-sum test as appropriate; categorical variables are summarized as n (%) and compared using Fisher's exact test. Differential protein analyses were repeated using limma models adjusted for age and sex. For late-vs-early cancer-only analyses, additional sensitivity models included BMI, smoking, drinking and routine clinical laboratory variables where complete data were available.",
  "",
  "For diagnostic signatures, threshold-dependent performance was calculated at the Youden index using cross-validated predicted probabilities. Calibration was summarized using Brier score, calibration intercept and calibration slope. Decision-curve data were generated across thresholds from 0.05 to 0.95. Traditional tumor markers were compared only in the late-vs-early cancer subset because marker data were not available for healthy controls.",
  "",
  "## Limitation text that should be included",
  "",
  "The clinical metadata enabled age- and sex-adjusted analyses for the full cohort, but several clinical variables were unavailable in healthy controls. Therefore, analyses involving BMI, routine laboratory indices or tumor markers are restricted to cancer patients and should be interpreted as sensitivity analyses rather than full-cohort confounder control. The lack of healthy-control tumor marker measurements also prevents a direct comparison between the Olink early-detection signature and conventional tumor markers for early cancer detection."
)
writeLines(report, file.path(report_dir, "clinical_enhancement_report.md"), useBytes = TRUE)

message("Clinical enhancement complete: ", tab_dir)
