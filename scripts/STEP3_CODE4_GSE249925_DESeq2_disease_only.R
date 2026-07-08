
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
# Detect sample column
# ----------------------------
possible_sample_cols <- c("gsm_id", "sample_id", "sample", "GSM", "geo_accession")
sample_col <- possible_sample_cols[possible_sample_cols %in% colnames(design)][1]

if (is.na(sample_col)) {
    stop("Could not detect sample ID column in design table.")
}

cat("\nDetected sample column: ", sample_col, "\n")

# ----------------------------
# Detect disease/group column
# ----------------------------
possible_group_cols <- c("disease_group", "group", "condition", "phenotype", "diagnosis", "hcm_status")
group_col <- possible_group_cols[possible_group_cols %in% colnames(design)][1]

if (is.na(group_col)) {
    stop("Could not detect disease/group column in design table.")
}

cat("Detected group column: ", group_col, "\n")

cat("\nRaw group counts:\n")
print(table(design[[group_col]], useNA="ifany"))

# ----------------------------
# Standardize group labels
# ----------------------------
design$group_for_de <- as.character(design[[group_col]])
design$group_for_de_lower <- tolower(design$group_for_de)

design$condition <- NA

design$condition[grepl("control|normal|non.failing|non-failing|donor", design$group_for_de_lower)] <- "Control"
design$condition[grepl("hcm|hypertrophic", design$group_for_de_lower)] <- "HCM"

if (any(is.na(design$condition))) {
    cat("\nRows with unresolved condition labels:\n")
    print(design[is.na(design$condition), c(sample_col, group_col)])
    stop("Some samples could not be assigned to HCM or Control.")
}

design$condition <- factor(design$condition, levels=c("Control", "HCM"))

cat("\nStandardized DE group counts:\n")
print(table(design$condition))

# ----------------------------
# Align samples
# ----------------------------
design_samples <- as.character(design[[sample_col]])
count_samples <- colnames(counts)

if (!setequal(design_samples, count_samples)) {
    cat("\nSamples in design but missing from count matrix:\n")
    print(setdiff(design_samples, count_samples))

    cat("\nSamples in count matrix but missing from design:\n")
    print(setdiff(count_samples, design_samples))

    stop("Design/count sample mismatch.")
}

counts <- counts[, design_samples]
stopifnot(all(colnames(counts) == as.character(design[[sample_col]])))

# ----------------------------
# Clean counts
# ----------------------------
counts <- as.matrix(counts)
mode(counts) <- "numeric"

# DESeq2 requires integer counts
counts[is.na(counts)] <- 0
counts[counts < 0] <- 0
counts <- round(counts)

# Remove genes with all-zero / very low counts
keep <- rowSums(counts) >= 10
counts_filtered <- counts[keep, ]

cat("\nGenes before count filter: ", nrow(counts), "\n")
cat("Genes after count filter: ", nrow(counts_filtered), "\n")

# ----------------------------
# Run DESeq2 disease-only model
# ----------------------------
dds <- DESeqDataSetFromMatrix(
    countData = counts_filtered,
    colData = design,
    design = ~ condition
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
res_df$model <- "DESeq2_disease_only"
res_df$n_hcm <- sum(design$condition == "HCM")
res_df$n_control <- sum(design$condition == "Control")

# Rename columns into project standard
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
# Save full DE table
# ----------------------------
full_out <- file.path(out_dir, "GSE249925_DESeq2_disease_only_DE_full_gene_level.csv")
write.csv(res_df, full_out, row.names=FALSE)

cat("\nFull DESeq2 gene-level table saved:\n")
cat(full_out, "\n")

# ----------------------------
# Extract locked MCI target genes
# ----------------------------
target_universe$gene_symbol_upper <- toupper(as.character(target_universe$gene_symbol))
res_df$gene_symbol_upper <- toupper(as.character(res_df$gene_symbol))

target_de <- target_universe %>%
    left_join(res_df, by="gene_symbol_upper", suffix=c("_target", "_de"))

if ("gene_symbol_target" %in% colnames(target_de)) {
    target_de$gene_symbol <- target_de$gene_symbol_target
}

target_out <- file.path(out_dir, "GSE249925_DESeq2_disease_only_DE_locked_MCI_target_genes.csv")
write.csv(target_de, target_out, row.names=FALSE)

cat("\nLocked MCI target-gene DE table saved:\n")
cat(target_out, "\n")

matched_count <- sum(!is.na(target_de$log2FC))
total_target_rows <- nrow(target_de)

cat("\nTarget gene matching:\n")
cat("Target universe disease-gene rows: ", total_target_rows, "\n")
cat("Rows matched in GSE249925 DE: ", matched_count, "\n")
cat("Rows not matched: ", total_target_rows - matched_count, "\n")

# ----------------------------
# Save normalized counts for QC/future plots
# ----------------------------
norm_counts <- counts(dds, normalized=TRUE)
norm_counts_out <- file.path(out_dir, "GSE249925_DESeq2_normalized_counts.csv")
write.csv(norm_counts, norm_counts_out)

# ----------------------------
# Save summary JSON manually
# ----------------------------
summary_lines <- c(
    "{",
    paste0('  "step": "MCI Step 3 Code 4 GSE249925 DESeq2 disease-only DE",'),
    paste0('  "timestamp": "', Sys.time(), '",'),
    paste0('  "cohort": "GSE249925",'),
    paste0('  "model": "DESeq2: counts ~ HCM_status",'),
    paste0('  "comparison": "HCM_vs_Control",'),
    paste0('  "log2FC_definition": "HCM minus Control",'),
    paste0('  "n_samples": ', nrow(design), ','),
    paste0('  "n_hcm": ', sum(design$condition == "HCM"), ','),
    paste0('  "n_control": ', sum(design$condition == "Control"), ','),
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
    paste0('  "normalized_counts_output": "', norm_counts_out, '"'),
    "}"
)

summary_path <- file.path(qc_dir, "GSE249925_DESeq2_disease_only_DE_summary.json")
writeLines(summary_lines, summary_path)

cat("\nSummary JSON saved:\n")
cat(summary_path, "\n")

# ----------------------------
# Print top results
# ----------------------------
cat("\nTop 25 DE genes by FDR:\n")
print(
    res_df %>%
        select(gene_symbol, log2FC, SE, test_stat, p_value, FDR, mean_normalized_count) %>%
        head(25)
)

cat("\nTop matched locked MCI target genes by FDR:\n")
target_display <- target_de %>%
    filter(!is.na(log2FC)) %>%
    arrange(FDR, p_value, desc(abs(log2FC)))

show_cols <- intersect(
    c("gene_symbol", "disease_group", "stratum", "is_prespecified_h4_gene",
      "present_in_strict_clinvar_filter", "log2FC", "SE", "p_value", "FDR"),
    colnames(target_display)
)

print(head(target_display[, show_cols], 40))

cat("\nCODE 4 COMPLETE — GSE249925 DESeq2 disease-only DE finished.\n")
