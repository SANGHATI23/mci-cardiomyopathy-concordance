
suppressPackageStartupMessages({
    library(limma)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(readr)
    library(dplyr)
    library(tibble)
})

expr_path <- "/content/drive/MyDrive/MCI_Project/data/processed/expression_matrices/GSE141910_expression_matrix_ensembl_values_HCM_Control.csv"
design_path <- "/content/drive/MyDrive/MCI_Project/data/processed/sample_metadata/GSE141910_analysis_design_table_available_samples.csv"
target_gene_path <- "/content/drive/MyDrive/MCI_Project/data/processed/clinvar/mci_target_gene_universe.csv"

out_dir <- "/content/drive/MyDrive/MCI_Project/results/differential_expression/GSE141910"
qc_dir <- "/content/drive/MyDrive/MCI_Project/results/quality_control/differential_expression"
log_dir <- "/content/drive/MyDrive/MCI_Project/results/logs"

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(qc_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(log_dir, recursive=TRUE, showWarnings=FALSE)

cat("Loading expression matrix...\n")
expr <- read.csv(expr_path, check.names=FALSE, row.names=1)

cat("Loading design table...\n")
design <- read.csv(design_path, check.names=FALSE)

cat("Loading target gene universe...\n")
target_universe <- read.csv(target_gene_path, check.names=FALSE)

cat("\nExpression dimensions:\n")
print(dim(expr))

cat("\nDesign dimensions:\n")
print(dim(design))

cat("\nDesign columns:\n")
print(colnames(design))

# ----------------------------
# 4. Validate design columns
# ----------------------------
required_cols <- c(
    "gsm_id",
    "disease_group",
    "disease_binary_HCM_vs_Control",
    "sex",
    "sex_binary_male_vs_female",
    "age_years"
)

missing_cols <- setdiff(required_cols, colnames(design))

if (length(missing_cols) > 0) {
    stop(paste("Missing required design columns:", paste(missing_cols, collapse=", ")))
}

# ----------------------------
# 5. Align expression columns to design rows
# ----------------------------
design$gsm_id <- as.character(design$gsm_id)

if (!setequal(design$gsm_id, colnames(expr))) {
    cat("\nSamples in design but missing from expression:\n")
    print(setdiff(design$gsm_id, colnames(expr)))

    cat("\nSamples in expression but missing from design:\n")
    print(setdiff(colnames(expr), design$gsm_id))

    stop("Expression/design sample mismatch.")
}

expr <- expr[, design$gsm_id]
stopifnot(all(colnames(expr) == design$gsm_id))

cat("\nSample alignment confirmed.\n")

# ----------------------------
# 6. Prepare phenotype and covariates
# ----------------------------
design$condition <- factor(design$disease_group, levels=c("Control", "HCM"))
design$age_years <- suppressWarnings(as.numeric(design$age_years))

design$sex_clean <- tolower(as.character(design$sex))
design$sex_clean[design$sex_clean %in% c("male", "m")] <- "Male"
design$sex_clean[design$sex_clean %in% c("female", "f")] <- "Female"
design$sex_clean[!(design$sex_clean %in% c("Male", "Female"))] <- NA
design$sex_clean <- factor(design$sex_clean)

complete_mask <- !is.na(design$condition) &
                 !is.na(design$age_years) &
                 !is.na(design$sex_clean)

design_complete <- design[complete_mask, ]
expr_complete <- expr[, design_complete$gsm_id]

cat("\nOriginal samples:", nrow(design), "\n")
cat("Complete samples for adjusted model:", nrow(design_complete), "\n")
cat("Excluded samples due to missing condition/age/sex:", nrow(design) - nrow(design_complete), "\n")

cat("\nCondition counts after filtering:\n")
print(table(design_complete$condition, useNA="ifany"))

cat("\nSex counts after filtering:\n")
print(table(design_complete$sex_clean, useNA="ifany"))

cat("\nAge summary after filtering:\n")
print(summary(design_complete$age_years))

if (sum(design_complete$condition == "HCM") < 2 || sum(design_complete$condition == "Control") < 2) {
    stop("Not enough HCM or Control samples for DE.")
}

if (length(unique(design_complete$sex_clean)) < 2) {
    stop("Sex has fewer than two levels after filtering.")
}

# ----------------------------
# 7. Clean expression matrix
# ----------------------------
expr_complete <- as.matrix(expr_complete)
mode(expr_complete) <- "numeric"

# Remove rows with any non-finite issue
row_ok <- rowSums(is.finite(expr_complete)) == ncol(expr_complete)
expr_complete <- expr_complete[row_ok, ]

# Remove zero-variance genes
row_sd <- apply(expr_complete, 1, sd)
expr_complete <- expr_complete[row_sd > 0, ]

cat("\nGenes after finite-value and variance filter:", nrow(expr_complete), "\n")

# ----------------------------
# 8. Ensembl to gene symbol mapping
# ----------------------------
ensembl_raw <- rownames(expr_complete)
ensembl_clean <- sub("\\..*$", "", ensembl_raw)

symbol_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys=unique(ensembl_clean),
    columns=c("SYMBOL"),
    keytype="ENSEMBL"
)

symbol_map <- symbol_map[!is.na(symbol_map$SYMBOL), ]
symbol_map <- symbol_map[!duplicated(symbol_map$ENSEMBL), ]

map_df <- data.frame(
    ensembl_gene_id = ensembl_raw,
    ensembl_clean = ensembl_clean,
    stringsAsFactors=FALSE
)

map_df <- map_df %>%
    left_join(symbol_map, by=c("ensembl_clean" = "ENSEMBL"))

map_df$gene_symbol <- ifelse(is.na(map_df$SYMBOL), map_df$ensembl_clean, map_df$SYMBOL)

cat("\nMapped gene symbols:", sum(!is.na(map_df$SYMBOL)), "of", nrow(map_df), "\n")

annotation_out <- file.path(out_dir, "GSE141910_ensembl_to_gene_symbol_annotation.csv")
write.csv(map_df, annotation_out, row.names=FALSE)

# ----------------------------
# 9. Run limma adjusted model
# ----------------------------
model_matrix <- model.matrix(
    ~ age_years + sex_clean + condition,
    data=design_complete
)

cat("\nModel matrix columns:\n")
print(colnames(model_matrix))

fit <- lmFit(expr_complete, model_matrix)
fit <- eBayes(fit)

coef_name <- "conditionHCM"

if (!(coef_name %in% colnames(model_matrix))) {
    stop(paste("Could not find coefficient:", coef_name))
}

coef_index <- which(colnames(model_matrix) == coef_name)

top <- topTable(
    fit,
    coef=coef_name,
    number=Inf,
    sort.by="P",
    adjust.method="BH"
)

top$ensembl_gene_id <- rownames(top)

# Standard errors for HCM coefficient
se <- fit$stdev.unscaled[, coef_index] * sqrt(fit$s2.post)
top$SE <- se[rownames(top)]

# Add project-standard columns
top <- top %>%
    left_join(map_df[, c("ensembl_gene_id", "gene_symbol")], by="ensembl_gene_id")

res_df <- top %>%
    transmute(
        feature_id = ensembl_gene_id,
        ensembl_gene_id = ensembl_gene_id,
        gene_symbol = gene_symbol,
        cohort = "GSE141910",
        disease = "HCM",
        comparison = "HCM_vs_Control",
        model = "limma_adjusted_age_sex",
        n_hcm = sum(design_complete$condition == "HCM"),
        n_control = sum(design_complete$condition == "Control"),
        n_samples_model = nrow(design_complete),
        mean_expression_hcm = rowMeans(expr_complete[, design_complete$condition == "HCM", drop=FALSE]),
        mean_expression_control = rowMeans(expr_complete[, design_complete$condition == "Control", drop=FALSE]),
        log2FC = logFC,
        SE = SE,
        moderated_t = t,
        p_value = P.Value,
        FDR = adj.P.Val,
        B_statistic = B,
        average_expression = AveExpr
    )

res_df$abs_log2FC <- abs(res_df$log2FC)

res_df <- res_df %>%
    arrange(FDR, p_value, desc(abs_log2FC))

# ----------------------------
# 10. Save full Ensembl-level DE
# ----------------------------
full_ensembl_out <- file.path(out_dir, "GSE141910_limma_adjusted_age_sex_DE_full_ensembl_level.csv")
write.csv(res_df, full_ensembl_out, row.names=FALSE)

cat("\nFull Ensembl-level adjusted DE table saved:\n")
cat(full_ensembl_out, "\n")

# ----------------------------
# 11. Create gene-symbol representative table
# ----------------------------
gene_level <- res_df %>%
    arrange(gene_symbol, FDR, p_value, desc(abs_log2FC)) %>%
    group_by(gene_symbol) %>%
    slice(1) %>%
    ungroup() %>%
    arrange(FDR, p_value, desc(abs_log2FC))

gene_level_out <- file.path(out_dir, "GSE141910_limma_adjusted_age_sex_DE_gene_symbol_level.csv")
write.csv(gene_level, gene_level_out, row.names=FALSE)

cat("\nGene-symbol-level adjusted DE table saved:\n")
cat(gene_level_out, "\n")

# ----------------------------
# 12. Extract locked MCI target genes
# ----------------------------
target_universe$gene_symbol_upper <- toupper(as.character(target_universe$gene_symbol))
gene_level$gene_symbol_upper <- toupper(as.character(gene_level$gene_symbol))

target_de <- target_universe %>%
    left_join(gene_level, by="gene_symbol_upper", suffix=c("_target", "_de"))

if ("gene_symbol_target" %in% colnames(target_de)) {
    target_de$gene_symbol <- target_de$gene_symbol_target
}

target_out <- file.path(out_dir, "GSE141910_limma_adjusted_age_sex_DE_locked_MCI_target_genes.csv")
write.csv(target_de, target_out, row.names=FALSE)

cat("\nLocked MCI target-gene adjusted DE table saved:\n")
cat(target_out, "\n")

matched_count <- sum(!is.na(target_de$log2FC))
total_target_rows <- nrow(target_de)

cat("\nTarget gene matching:\n")
cat("Target universe disease-gene rows:", total_target_rows, "\n")
cat("Rows matched in GSE141910 adjusted DE:", matched_count, "\n")
cat("Rows not matched:", total_target_rows - matched_count, "\n")

# ----------------------------
# 13. Save final MCI DE copies
# ----------------------------
final_full_out <- file.path(out_dir, "GSE141910_FINAL_DE_for_MCI_gene_symbol_level.csv")
final_target_out <- file.path(out_dir, "GSE141910_FINAL_DE_for_MCI_locked_MCI_target_genes.csv")

write.csv(gene_level, final_full_out, row.names=FALSE)
write.csv(target_de, final_target_out, row.names=FALSE)

# ----------------------------
# 14. Save summary JSON
# ----------------------------
summary_lines <- c(
    "{",
    paste0('  "step": "MCI Step 3 Code 10 GSE141910 limma adjusted DE",'),
    paste0('  "timestamp": "', Sys.time(), '",'),
    paste0('  "cohort": "GSE141910",'),
    paste0('  "final_model_for_MCI": "limma adjusted age sex",'),
    paste0('  "final_model_formula": "expression ~ age_years + sex + condition",'),
    paste0('  "comparison": "HCM_vs_Control",'),
    paste0('  "log2FC_definition": "Adjusted HCM minus Control",'),
    paste0('  "n_samples_original": ', nrow(design), ','),
    paste0('  "n_samples_adjusted_model": ', nrow(design_complete), ','),
    paste0('  "n_samples_excluded_missing_covariates": ', nrow(design) - nrow(design_complete), ','),
    paste0('  "n_hcm": ', sum(design_complete$condition == "HCM"), ','),
    paste0('  "n_control": ', sum(design_complete$condition == "Control"), ','),
    paste0('  "n_ensembl_features_tested": ', nrow(res_df), ','),
    paste0('  "n_gene_symbol_rows": ', nrow(gene_level), ','),
    paste0('  "n_fdr_lt_0_05": ', sum(gene_level$FDR < 0.05, na.rm=TRUE), ','),
    paste0('  "n_fdr_lt_0_10": ', sum(gene_level$FDR < 0.10, na.rm=TRUE), ','),
    paste0('  "n_nominal_p_lt_0_05": ', sum(gene_level$p_value < 0.05, na.rm=TRUE), ','),
    paste0('  "target_universe_rows": ', total_target_rows, ','),
    paste0('  "target_rows_matched_in_de": ', matched_count, ','),
    paste0('  "target_rows_not_matched_in_de": ', total_target_rows - matched_count, ','),
    paste0('  "full_ensembl_de_output": "', full_ensembl_out, '",'),
    paste0('  "gene_symbol_de_output": "', gene_level_out, '",'),
    paste0('  "target_gene_de_output": "', target_out, '",'),
    paste0('  "final_full_de_output": "', final_full_out, '",'),
    paste0('  "final_target_gene_de_output": "', final_target_out, '",'),
    paste0('  "annotation_output": "', annotation_out, '"'),
    "}"
)

summary_path <- file.path(qc_dir, "GSE141910_FINAL_DE_for_MCI_summary.json")
writeLines(summary_lines, summary_path)

cat("\nSummary JSON saved:\n")
cat(summary_path, "\n")

# ----------------------------
# 15. Print top results
# ----------------------------
cat("\nTop 25 adjusted DE genes by FDR:\n")
print(
    gene_level %>%
        select(gene_symbol, ensembl_gene_id, log2FC, SE, moderated_t, p_value, FDR, average_expression) %>%
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

cat("\nCODE 10 COMPLETE — GSE141910 limma adjusted DE finished and locked for MCI.\n")
