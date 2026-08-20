# Supplementary Code 1

This archive contains the analysis scripts used for the Olink Reveal gastric cancer plasma proteomics manuscript.

## Contents

- `scripts/01_olink_primary_analysis.R`: primary Olink quality-control, differential, trend, enrichment and modelling analyses.
- `scripts/02_tcga_stad_external_validation.R`: TCGA-STAD tissue transcriptomic support analysis using public UCSC Xena/GDC files.
- `scripts/03_export_individual_panels.R`: individual figure-panel export.
- `scripts/04_clinical_enhancement_analysis.R`: clinical metadata integration, baseline tables and covariate-adjusted sensitivity analyses.
- `scripts/05_export_clinical_review_panels.R`: clinical robustness figure-panel export.
- `scripts/06_insert_figures_into_manuscript_docx.py`: insertion of composite figures into the manuscript DOCX.
- `scripts/07_finalize_manuscript_tables_declarations.py`: final manuscript table and declaration updates.

## Data availability

Raw participant-level Olink NPX and clinical metadata are not included in this code archive because they are subject to institutional and participant-privacy restrictions. The scripts are provided to document the analysis workflow and can be rerun after updating local input and output paths.

The TCGA-STAD public data used by `scripts/02_tcga_stad_external_validation.R` are downloaded from the UCSC Xena GDC hub:

- STAR TPM expression matrix: https://gdc.xenahubs.net/download/TCGA-STAD.star_tpm.tsv.gz
- Clinical annotations: https://gdc.xenahubs.net/download/TCGA-STAD.clinical.tsv.gz
- Survival annotations: https://gdc.xenahubs.net/download/TCGA-STAD.survival.tsv.gz

## Software

The analyses were performed in R using packages including `readxl`, `dplyr`, `tidyr`, `limma`, `glmnet`, `pROC`, `clusterProfiler`, `org.Hs.eg.db`, `ggplot2`, `ggrepel`, `patchwork`, `svglite`, `ragg`, `data.table`, `curl` and `survival`.

## Reuse note

Before reuse, set the following environment variables or edit the corresponding defaults in the scripts:

- `GC_OLINK_SUMMARY_DIR`: directory containing the Olink output summary files.
- `GC_OLINK_CLINICAL_FILE`: clinical metadata workbook.
- `GC_OLINK_RESULTS_DIR`: output directory for tables, figures, reports and cached public data.
- `GC_OLINK_PROJECT_DIR`: project directory used by figure/document export helper scripts.
