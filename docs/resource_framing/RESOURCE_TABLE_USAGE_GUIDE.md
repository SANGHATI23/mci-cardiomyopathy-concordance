
# Resource Table Usage Guide

## Master table

File:

results/resource_tables/MCI_CARDIOMYOPATHY_RESOURCE_MASTER_TABLE_v0_2.csv

This is the main browser-ready table for the MCI cardiomyopathy resource.

## Recommended user workflow

### 1. Search by gene

Use gene_symbol to find a ClinVar-annotated cardiomyopathy gene.

### 2. Check primary MCI

Use:

- MCI
- MCI_tier
- MCI_CI95_lower
- MCI_CI95_upper
- bootstrap_majority_tier

These fields show whether the gene has reproducible transcript-level evidence across current HCM cohorts.

### 3. Check GTEx baseline status

Use:

- sigma_disease
- sigma_GTEx
- sigma_disease_to_GTEx_ratio
- GTEx_low_confidence_flag
- GTEx_adjustment_status

These fields show whether disease-associated variability exceeds normal left-ventricle expression variability.

### 4. Check DCM generalization

Use:

- has_DCM_generalization_entry
- DCM_MCI_if_available
- DCM_tier_if_available

These fields show whether a comparable DCM evidence layer was available.

### 5. Read plain-language interpretation

Use:

- resource_interpretation

This field summarizes how to interpret the row for target-evidence review.

## Important interpretation examples

### HIGH MCI plus GTEx-low-confidence

This means the gene is reproducible across HCM cohorts, but the magnitude of disease-associated variability does not exceed GTEx normal left-ventricle baseline variability. The gene should be treated as transcript-concordant but not as an unqualified transcript-level target.

### UNSTABLE MCI

This means the gene has inconsistent transcript-level evidence across current HCM cohorts. It may still be clinically important, especially if its mechanism is protein-level, structural, splicing-level, or cell-type-specific.

### DCM MCI missing

This does not mean the gene has no DCM relevance. It means DCM generalization evidence was not available or not attached in the current resource version.

## Suggested manuscript wording

The MCI resource should be described as an evidence-auditing layer that distinguishes clinical pathogenicity from transcriptomic reproducibility and normal-tissue variability. It is not a clinical decision tool and should not be interpreted as evidence of therapeutic validity without additional validation.
