# STEP 10 — Create README.md for GitHub repository

from pathlib import Path

PROJECT_ROOT = Path("/content/drive/MyDrive/MCI_Project")

readme_text = """# Molecular Concordance Index Cardiomyopathy Project

## Project Overview

This repository implements the Molecular Concordance Index (MCI) framework for evaluating whether ClinVar-annotated cardiomyopathy genes show reproducible phenotype-associated transcriptomic signatures across independent public datasets.

The project focuses on hypertrophic cardiomyopathy (HCM) and dilated cardiomyopathy (DCM), using public GEO transcriptomic datasets, ClinVar gene annotations, and GTEx left ventricle baseline expression variability.

## Core Research Question

Do clinically annotated pathogenic or likely pathogenic cardiomyopathy genes demonstrate stable and reproducible transcriptomic behavior across independent cohorts?

## Main Hypothesis

Sarcomeric genes such as MYH7, MYBPC3, TNNT2, and related genes are expected to show lower transcript-level concordance than non-sarcomeric or RNA-processing genes such as LMNA, RBM20, PLN, FLNC, and BAG3.

## Molecular Concordance Index Formula

MCI_g = 0.40 * D_g + 0.35 * S_g + 0.25 * R_g

Where:

- D_g = direction agreement across cohorts
- S_g = effect size consistency across cohorts
- R_g = statistical reproducibility across cohorts

## MCI Tier Thresholds

| Tier | Threshold | Interpretation |
|---|---|---|
| High | MCI >= 0.70 | Stable phenotype-associated transcriptomic signature |
| Moderate | 0.45 <= MCI < 0.70 | Partial or conditional concordance |
| Unstable | MCI < 0.45 | Inconsistent or cohort-dependent expression pattern |

## Current Repository Status

This repository currently contains the initial executable project scaffold:

- Project configuration file
- Dataset manifest
- Dataset eligibility validation
- Pre-specified target gene manifest
- Reusable MCI scoring script
- Mock differential-expression input
- Mock MCI output table

## Planned Pipeline

1. Download and curate public GEO datasets
2. Extract and filter ClinVar cardiomyopathy P/LP genes
3. Harmonize gene identifiers
4. Run per-cohort differential-expression analysis
5. Compute MCI scores
6. Add GTEx baseline variability adjustment
7. Run sensitivity analyses
8. Build interactive concordance browser
9. Deposit final results to GitHub and Zenodo

## Repository Structure

```text
MCI_Project/
├── data/
│   ├── raw/
│   ├── processed/
│   └── metadata/
├── results/
│   ├── differential_expression/
│   ├── mci_scores/
│   └── figures/
├── scripts/
│   └── mci_scoring.py
├── logs/
├── project_config.json
└── README.md
