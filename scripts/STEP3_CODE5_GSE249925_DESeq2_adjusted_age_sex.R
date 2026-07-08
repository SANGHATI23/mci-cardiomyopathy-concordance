
suppressPackageStartupMessages({
    library(DESeq2)
    library(readr)
    library(dplyr)
    library(tibble)
})

count_path <- "/content/drive/MyDrive/MCI_Project/data/processed/expression_matrices/GSE249925_count_matrix_gene_symbol_raw_counts.csv"
design_path <- "/content/drive/MyDrive/MCI_Project/data/processed/sample_metadata/GSE249925_analysis_design_table.csv"
target_gene_path <- "/content/drive/MyDrive/MCI_Project/data/processed/clinvar/mci_target_gene_universe.csv"

out_dir <- "/content/drive/MyDrive/MCI_Project/results/differential_expression/GSE249925"
qc_dir <- "/content/drive/MyDrive/MCI_Project/results/quality_control/differential_expression"
log_dir <- "/content/drive/MyDrive/MCI_Project/results/logs"

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(qc_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(log_dir, recursive=TRUE, showWarnings=FALSE)

cat("Loading GSE249925 count matrix...\n")
counts <- read.csv(count_path, check.names=FALSE, row.names=1)

cat("Loading design table...\n")
design <- read.csv(design_path, check.names=FALSE)

cat("Loading target gene universe...\n")
target_universe <- read.csv(target_gene_path, check.names=FALSE)

cat("\nCount matrix dimensions:\n")
print(dim(counts))

cat("\nDesign dimensions:\n")
print(dim(design))

cat("\nDesign columns:\n")
print(colnames(design))

# ----------------------------
# 3. Required design columns
# ----------------------------
required_cols <- c(
    "gsm_id",
    "disease_group",
    "age_years",
    "sex",
    "sex_binary_male_vs_female"
)

missing_cols <- setdiff(required_cols, colnames(design))

if (length(missing_cols) > 0) {
    stop(paste("Missing required design columns:", paste(missing_cols, collapse=", ")))
}

# ----------------------------
# 4. Standardize disease group
# ----------------------------
design$group_for_de <- as.character(design$disease_group)
design$group_for_de_lower <- tolower(design$group_for_de)

design$condition <- NA
design$condition[grepl("control|normal|non.failing|non-failing|donor", design$group_for_de_lower)] <- "Control"
design$condition[grepl("hcm|hypertrophic", design$group_for_de_lower)] <- "HCM"

if (any(is.na(design$condition))) {
    cat("\nUnresolved condition rows:\n")
    print(design[is.na(design$condition), c("gsm_id", "disease_group")])
    stop("Some samples could not be assigned to HCM or Control.")
}

design$condition <- factor(design$condition, levels=c("Control", "HCM"))

cat("\nCondition counts before covariate filtering:\n")
print(table(design$condition, useNA="ifany"))

# ----------------------------
# 5. Clean covariates
# ----------------------------
design$age_years <- suppressWarnings(as.numeric(design$age_years))

# Use sex text column if available; keep clean factor
design$sex_clean <- tolower(as.character(design$sex))
design$sex_clean[design$sex_clean %in% c("male", "m")] <- "Male"
design$sex_clean[design$sex_clean %in% c("female", "f")] <- "Female"
design$sex_clean[!(design$sex_clean %in% c("Male", "Female"))] <- NA
design$sex_clean <- factor(design$sex_clean)

cat("\nAge missing count:\n")
print(sum(is.na(design$age_years)))

cat("\nSex missing count:\n")
print(sum(is.na(design$sex_clean)))

cat("\nSex counts before covariate filtering:\n")
print(table(design$sex_clean, useNA="ifany"))

# ----------------------------
# 6. Keep complete covariate samples only
# ----------------------------
complete_mask <- !is.na(design$gsm_id) &
                 !is.na(design$condition) &
                 !is.na(design$age_years) &
                 !is.na(design$sex_clean)

design_complete <- design[complete_mask, ]

cat("\nOriginal samples:", nrow(design), "\n")
cat("Complete samples for adjusted model:", nrow(design_complete), "\n")
cat("Excluded samples due to missing condition/age/sex:", nrow(design) - nrow(design_complete), "\n")

cat("\nCondition counts after covariate filtering:\n")
print(table(design_complete$condition, useNA="ifany"))

cat("\nSex counts after covariate filtering:\n")
print(table(design_complete$sex_clean, useNA="ifany"))

cat("\nAge summary after covariate filtering:\n")
print(summary(design_complete$age_years))

if (sum(design_complete$condition == "HCM") < 2 || sum(design_complete$condition == "Control") < 2) {
    stop("Not enough HCM or Control samples after covariate filtering.")
}

if (length(unique(design_complete$sex_clean)) < 2) {
    stop("Sex has fewer than two levels after filtering; adjusted model cannot include sex.")
}

# ----------------------------
# 7. Align counts to complete design
# ----------------------------
design_samples <- as.character(design_complete$gsm_id)
count_samples <- colnames(counts)

if (!setequal(design_samples, count_samples)) {
    cat("\nSamples in complete design but missing from count matrix:\n")
    print(setdiff(design_samples, count_samples))

    cat("\nSamples in count matrix but not used in complete adjusted model:\n")
    print(setdiff(count_samples, design_samples))
}

counts <- counts[, design_samples]
stopifnot(all(colnames(counts) == as.character(design_complete$gsm_id)))

# ----------------------------
# 8. Clean counts
# ----------------------------
counts <- as.matrix(counts)
mode(counts) <- "numeric"

counts[is.na(counts)] <- 0
counts[counts < 0] <- 0
counts <- round(counts)

keep <- rowSums(counts) >= 10
counts_filtered <- counts[keep, ]

cat("\nGenes before count filter:", nrow(counts), "\n")
cat("Genes after count filter:", nrow(counts_filtered), "\n")

# ----------------------------
# 9. Run adjusted DESeq2
# ----------------------------
dds <- DESeqDataSetFromMatrix(
    countData = counts_filtered,
    colData = design_complete,
    design = ~ age_years + sex_clean + condition
)

dds <- DESeq(dds)

res <- results(
    dds,
    contrast=c("condition", "HCM", "Control"),
    alpha=0.05
)

res_df <- as.data.frame(res)
res_df$gene_symbol <- rownames(res_df)
res_df$feature_id <- rownames(res_df)
res_df$cohort <- "GSE249925"
res_df$disease <- "HCM"
res_df$comparison <- "HCM_vs_Control"
res_df$model <- "DESeq2_adjusted_age_sex"
res_df$n_hcm <- sum(design_complete$condition == "HCM")
res_df$n_control <- sum(design_complete$condition == "Control")
res_df$n_samples_model <- nrow(design_complete)

res_df <- res_df %>%
    rename(
        mean_normalized_count = baseMean,
        log2FC = log2FoldChange,
        SE = lfcSE,
        test_stat = stat,
        p_value = pvalue,
        FDR = padj
    )

res_df$abs_log2FC <- abs(res_df$log2FC)

res_df <- res_df %>%
    arrange(FDR, p_value, desc(abs_log2FC))

# ----------------------------
# 10. Save full adjusted DE table
# ----------------------------
full_out <- file.path(out_dir, "GSE249925_DESeq2_adjusted_age_sex_DE_full_gene_level.csv")
write.csv(res_df, full_out, row.names=FALSE)

cat("\nFull adjusted DESeq2 table saved:\n")
cat(full_out, "\n")

# ----------------------------
# 11. Extract locked MCI target genes
# ----------------------------
target_universe$gene_symbol_upper <- toupper(as.character(target_universe$gene_symbol))
res_df$gene_symbol_upper <- toupper(as.character(res_df$gene_symbol))

target_de <- target_universe %>%
    left_join(res_df, by="gene_symbol_upper", suffix=c("_target", "_de"))

if ("gene_symbol_target" %in% colnames(target_de)) {
    target_de$gene_symbol <- target_de$gene_symbol_target
}

target_out <- file.path(out_dir, "GSE249925_DESeq2_adjusted_age_sex_DE_locked_MCI_target_genes.csv")
write.csv(target_de, target_out, row.names=FALSE)

cat("\nLocked MCI target-gene adjusted DE table saved:\n")
cat(target_out, "\n")

matched_count <- sum(!is.na(target_de$log2FC))
total_target_rows <- nrow(target_de)

cat("\nTarget gene matching:\n")
cat("Target universe disease-gene rows:", total_target_rows, "\n")
cat("Rows matched in GSE249925 adjusted DE:", matched_count, "\n")
cat("Rows not matched:", total_target_rows - matched_count, "\n")

# ----------------------------
# 12. Save normalized counts
# ----------------------------
norm_counts <- counts(dds, normalized=TRUE)
norm_counts_out <- file.path(out_dir, "GSE249925_DESeq2_adjusted_age_sex_normalized_counts.csv")
write.csv(norm_counts, norm_counts_out)

# ----------------------------
# 13. Save excluded sample list
# ----------------------------
excluded_samples <- design[!complete_mask, ]
excluded_out <- file.path(qc_dir, "GSE249925_adjusted_age_sex_excluded_samples.csv")
write.csv(excluded_samples, excluded_out, row.names=FALSE)

# ----------------------------
# 14. Save summary JSON
# ----------------------------
summary_lines <- c(
    "{",
    paste0('  "step": "MCI Step 3 Code 5 GSE249925 DESeq2 covariate-adjusted DE",'),
    paste0('  "timestamp": "', Sys.time(), '",'),
    paste0('  "cohort": "GSE249925",'),
    paste0('  "model": "DESeq2: counts ~ age_years + sex + HCM_status",'),
    paste0('  "comparison": "HCM_vs_Control",'),
    paste0('  "log2FC_definition": "Adjusted HCM minus Control",'),
    paste0('  "n_samples_original": ', nrow(design), ','),
    paste0('  "n_samples_adjusted_model": ', nrow(design_complete), ','),
    paste0('  "n_samples_excluded_missing_covariates": ', nrow(design) - nrow(design_complete), ','),
    paste0('  "n_hcm": ', sum(design_complete$condition == "HCM"), ','),
    paste0('  "n_control": ', sum(design_complete$condition == "Control"), ','),
    paste0('  "n_genes_before_filter": ', nrow(counts), ','),
    paste0('  "n_genes_after_filter": ', nrow(counts_filtered), ','),
    paste0('  "n_fdr_lt_0_05": ', sum(res_df$FDR < 0.05, na.rm=TRUE), ','),
    paste0('  "n_fdr_lt_0_10": ', sum(res_df$FDR < 0.10, na.rm=TRUE), ','),
    paste0('  "n_nominal_p_lt_0_05": ', sum(res_df$p_value < 0.05, na.rm=TRUE), ','),
    paste0('  "target_universe_rows": ', total_target_rows, ','),
    paste0('  "target_rows_matched_in_de": ', matched_count, ','),
    paste0('  "target_rows_not_matched_in_de": ', total_target_rows - matched_count, ','),
    paste0('  "full_de_output": "', full_out, '",'),
    paste0('  "target_gene_de_output": "', target_out, '",'),
    paste0('  "normalized_counts_output": "', norm_counts_out, '",'),
    paste0('  "excluded_samples_output": "', excluded_out, '"'),
    "}"
)

summary_path <- file.path(qc_dir, "GSE249925_DESeq2_adjusted_age_sex_DE_summary.json")
writeLines(summary_lines, summary_path)

cat("\nSummary JSON saved:\n")
cat(summary_path, "\n")

# ----------------------------
# 15. Print top results
# ----------------------------
cat("\nTop 25 adjusted DE genes by FDR:\n")
print(
    res_df %>%
        select(gene_symbol, log2FC, SE, test_stat, p_value, FDR, mean_normalized_count) %>%
        head(25)
)

cat("\nTop matched locked MCI target genes by adjusted FDR:\n")
target_display <- target_de %>%
    filter(!is.na(log2FC)) %>%
    arrange(FDR, p_value, desc(abs(log2FC)))

show_cols <- intersect(
    c("gene_symbol", "disease_group", "stratum", "is_prespecified_h4_gene",
      "present_in_strict_clinvar_filter", "log2FC", "SE", "p_value", "FDR"),
    colnames(target_display)
)

print(head(target_display[, show_cols], 40))

cat("\nCODE 5 COMPLETE — GSE249925 covariate-adjusted DESeq2 finished.\n")
