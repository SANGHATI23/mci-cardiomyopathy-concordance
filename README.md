# Molecular Concordance Index Cardiomyopathy Project

Reproducible computational pipeline for Molecular Concordance Index scoring of ClinVar-annotated cardiomyopathy genes across public HCM/DCM transcriptomic cohorts.

## Project Goal

This project evaluates whether clinically annotated cardiomyopathy genes show stable phenotype-associated transcriptomic signatures across independent cohorts.


MCI evaluates whether clinically annotated cardiomyopathy genes
show reproducible transcriptomic evidence across independent cohorts.

## Start Here

For a step-by-step reproduction and demonstration of the project:

1. [Project Setup and Data Acquisition](notebooks/01_MCI_Project_Setup_and_Data_Acquisition.ipynb)
2. [GTEx Baseline Adjustment](notebooks/02_GTEx_Baseline_Adjustment.ipynb)
3. [DCM Generalization and External Validation](notebooks/03_DCM_Generalization_and_External_Validation.ipynb)
4. [MCI Resource and Browser](notebooks/04_MCI_Resource_and_Browser.ipynb)
5. [Benchmarking and Sensitivity Analysis](notebooks/05_MCI_Benchmarking_and_Sensitivity_Analysis.ipynb)

   Clinical gene evidence
        ↓
Eligible HCM cohorts
        ↓
Cohort-specific differential expression
        ↓
Cross-cohort concordance
        ↓
MCI
        ↓
GTEx baseline adjustment
        ↓
Bootstrap uncertainty
        ↓
DCM portability
        ↓
Held-out/GWAS validation
        ↓
Benchmarking & sensitivity
        ↓
Browsable MCI resource

## MCI Formula

`MCI_g = 0.40 * D_g + 0.35 * S_g + 0.25 * R_g`

- D_g = direction agreement
- S_g = effect size consistency
- R_g = statistical reproducibility

## Current Status

- Project configuration created
- Dataset manifest created
- Dataset eligibility validation created
- Target gene manifest created
- Reusable MCI scoring script created
- Mock DE and MCI output generated

## Author

Sanghati Basu


<!-- TASK11_START -->
## Task 11: Benchmarking and Sensitivity Analysis

The MCI resource now includes formal benchmarking against equal-weight
scoring, direction voting, Fisher combined significance, signed Stouffer
combined significance, and DerSimonian-Laird random-effects
meta-analysis.

Task 11 also includes:

- alternative-weight and leave-one-component-out sensitivity analyses;
- MODERATE and HIGH threshold-shift analyses;
- leave-one-HCM-cohort-out recomputation;
- held-out GSE160997 expression benchmarking;
- GWAS Catalog convergence benchmarking;
- publication-ready tables and vector figures.

### Task 11 resources

- [Reproducibility instructions](results/task11_benchmarking/README.md)
- [Methods and results summary](results/task11_benchmarking/TASK11_METHODS_AND_RESULTS_SUMMARY.md)
- [Input-data manifest with SHA-256 checksums](results/task11_benchmarking/TASK11_INPUT_DATA_MANIFEST.csv)
- [Google Colab notebook](notebooks/TASK11_MCI_BENCHMARKING_AND_SENSITIVITY_ANALYSIS.ipynb)
- [Publication-ready outputs](results/task11_benchmarking/publication_ready/)
- [Publication-output manifest](results/task11_benchmarking/publication_ready/TASK11_PUBLICATION_OUTPUT_MANIFEST.csv)

### Task 11 environment

- `requirements-task11.txt`
- `environment-task11.yml`

### Interpretation boundary

The benchmarking supports MCI as a transparent multidimensional
transcriptomic evidence-audit framework. It does not establish MCI as a
clinically validated predictor, demonstrate universal superiority over
all simpler methods, or show complete cohort-invariant tier assignment.
<!-- TASK11_END -->
