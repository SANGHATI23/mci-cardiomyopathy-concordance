# Molecular Concordance Index Cardiomyopathy Project

Reproducible computational pipeline for Molecular Concordance Index scoring of ClinVar-annotated cardiomyopathy genes across public HCM/DCM transcriptomic cohorts.

## Project Goal

This project evaluates whether clinically annotated cardiomyopathy genes show stable phenotype-associated transcriptomic signatures across independent cohorts.

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
