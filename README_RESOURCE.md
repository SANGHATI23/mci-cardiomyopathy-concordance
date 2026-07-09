
# Molecular Concordance Index Cardiomyopathy Resource

## Overview

This repository contains the Molecular Concordance Index, or MCI, cardiomyopathy resource. It is a reproducible transcriptomic evidence-auditing framework for ClinVar-annotated hypertrophic cardiomyopathy and dilated cardiomyopathy genes.

The goal of this resource is not to claim that every clinically annotated cardiomyopathy gene has reproducible bulk-transcriptomic evidence. Instead, the resource makes that question explicit and queryable.

MCI separates four evidence layers that are often conflated:

1. Clinical pathogenicity annotation.
2. Cross-cohort transcriptomic reproducibility.
3. Normal left-ventricle expression variability.
4. Independent validation or convergence evidence.

The central resource table allows users to inspect whether a gene is transcript-concordant, transcript-unstable, GTEx-low-confidence, DCM-generalizable, or insufficiently covered by current public data.

## Current resource version

Version: v0.2-resource-reframe-corrected-source

Primary master table:

results/resource_tables/MCI_CARDIOMYOPATHY_RESOURCE_MASTER_TABLE_v0_2.csv

Data dictionary:

results/resource_tables/MCI_CARDIOMYOPATHY_RESOURCE_DATA_DICTIONARY_v0_2.csv

## Current resource content

The current HCM-primary resource contains:

- 49 ClinVar-annotated HCM genes
- 18 HIGH MCI genes
- 8 MODERATE MCI genes
- 21 UNSTABLE MCI genes
- 2 genes with insufficient cohort coverage
- 47 genes with disease-to-GTEx variability ratios
- 46 genes flagged as GTEx-low-confidence
- 3 genes not flagged as GTEx-low-confidence
- 22 genes with attached DCM generalization entries

## How to interpret MCI

MCI is a composite transcriptomic concordance score:

MCI_g = 0.40 * D_g + 0.35 * S_g + 0.25 * R_g

where:

- D_g = direction agreement across cohorts
- S_g = effect-size consistency across cohorts
- R_g = statistical reproducibility across cohorts

Tier labels:

- HIGH: MCI >= 0.70
- MODERATE: 0.45 <= MCI < 0.70
- UNSTABLE: MCI < 0.45
- INSUFFICIENT_COHORT_COVERAGE: fewer than two eligible cohorts

## GTEx baseline interpretation

The GTEx layer compares disease-associated variability against normal human left-ventricle expression variability.

ratio = sigma_disease / sigma_GTEx

If the ratio is less than or equal to 1.0, the gene is flagged as GTEx-low-confidence. This does not mean the gene is biologically irrelevant. It means the disease-associated bulk-transcriptomic variability does not exceed normal GTEx left-ventricle baseline variability in the current data.

This distinction is important because a gene may be clinically pathogenic and transcript-concordant while still requiring orthogonal support from genotype-stratified cohorts, single-cell analysis, proteomics, or functional validation.

## Main claim boundary

This resource does not claim that MCI proves causality, clinical actionability, or therapeutic validity.

The resource is designed for evidence auditing and target-prioritization support. It should be used to identify genes whose transcript-level evidence is reproducible, unstable, baseline-confounded, disease-context-specific, or underpowered.

## Manuscript framing

This project is best framed as a database/resource paper rather than as a high-impact biological discovery paper.

Recommended manuscript title:

Molecular Concordance Index: a reproducible transcriptomic resource for evaluating ClinVar-annotated cardiomyopathy genes across human HCM and DCM cohorts

## Repository structure

results/resource_tables/
  MCI_CARDIOMYOPATHY_RESOURCE_MASTER_TABLE_v0_2.csv
  MCI_CARDIOMYOPATHY_RESOURCE_DATA_DICTIONARY_v0_2.csv

manuscript/resource_paper/
  00_RESOURCE_PAPER_POSITIONING.md
  01_REVISED_RESOURCE_ABSTRACT.md
  02_REVISED_CONTRIBUTIONS.md
  03_RESOURCE_MANUSCRIPT_OUTLINE.md

docs/resource_framing/
  CLAIM_BOUNDARY_STATEMENT.md
  RESOURCE_TABLE_USAGE_GUIDE.md
  NON_SIGNIFICANT_RESULTS_FRAMING.md

scripts/scrna_contextualization/
  Placeholder for future public single-cell RNA-seq contextualization module.

shiny_app/
  Browser files or future browser outputs.

## Suggested citation language

This repository provides a reproducible MCI resource for evaluating transcriptomic concordance of ClinVar-annotated cardiomyopathy genes across public human heart datasets. Users should cite the associated manuscript or preprint once available.

## Status

Active resource-paper reframing in progress.
