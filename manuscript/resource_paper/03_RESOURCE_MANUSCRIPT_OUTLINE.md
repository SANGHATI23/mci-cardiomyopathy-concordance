# Resource-Paper Manuscript Outline

## Title

Molecular Concordance Index: a reproducible transcriptomic resource for evaluating ClinVar-annotated cardiomyopathy genes across human HCM and DCM cohorts

## Abstract

Use the revised resource-paper abstract in `01_REVISED_RESOURCE_ABSTRACT.md`.

## Introduction

1. Cardiomyopathy genetics depends heavily on ClinVar and disease-gene annotation.
2. Clinical pathogenicity does not guarantee reproducible transcriptomic perturbation.
3. Public bulk transcriptomic cohorts provide an opportunity to audit transcript-level reproducibility.
4. Current gaps:
   - no reusable gene-level concordance resource for ClinVar cardiomyopathy genes,
   - limited separation of disease-cohort reproducibility from normal-tissue variability,
   - limited transparent reporting of underpowered or non-significant validation attempts,
   - limited browser-ready tools for gene-level evidence review.
5. This study presents MCI as a reproducible resource.

## Results

### 1. Resource construction and cohort coverage
Report HCM and DCM datasets, ClinVar gene universe, and score-eligible genes.

### 2. Primary HCM MCI resource
Report 49 HCM genes scored, tier counts, and top high/moderate/unstable genes.

### 3. Bootstrap uncertainty and tier stability
Report 1,000-iteration bootstrap, 47 genes with CIs, wide intervals, and tier-majority behavior.

### 4. GTEx baseline benchmarking
Report GTEx v8 Heart - Left Ventricle sample extraction, 48 of 49 GTEx matches, and low-confidence flags.

### 5. Mechanism-stratified use case
Report sarcomeric vs non-sarcomeric comparison as a pre-specified demonstration analysis, not as the main discovery claim.

### 6. DCM generalization module
Report DCM results as disease-context evaluation. Emphasize that the expected pattern was not confirmed.

### 7. Held-out and GWAS validation modules
Report directionally supportive but non-significant held-out and GWAS results.

### 8. Browser-ready resource outputs
Describe downloadable tables, figures, and planned or implemented Shiny browser.

### 9. Optional scRNA-seq contextualization module
If completed, add cell-type enrichment/expression context for MCI genes.

## Methods

1. Dataset acquisition and eligibility
2. ClinVar HCM/DCM gene universe construction
3. Differential-expression processing
4. Harmonization and batch-aware processing
5. MCI formula
6. Tier assignment
7. Bootstrap confidence intervals
8. GTEx baseline adjustment
9. DCM generalization
10. Held-out validation
11. GWAS convergence
12. Shiny/browser implementation
13. Optional scRNA-seq contextualization

## Discussion

1. MCI as a resource for translational evidence auditing.
2. Why clinical pathogenicity and transcript reproducibility are non-equivalent.
3. How GTEx benchmarking changes interpretation of high-MCI genes.
4. What non-significant validation modules mean in a resource-paper context.
5. Use cases for target selection, biomarker review, and study design.
6. Limitations:
   - few cohorts,
   - bulk tissue only,
   - missing genotype resolution,
   - no wet-lab validation,
   - incomplete scRNA/proteomic context.
7. Future directions:
   - scRNA-seq,
   - genotype-stratified MCI,
   - proteomics,
   - larger biobank-scale validation,
   - browser expansion.

## Data and Code Availability

Include GitHub, Zenodo, Shiny URL, and versioned CSV table locations once finalized.
