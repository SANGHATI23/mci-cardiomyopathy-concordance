# MCI Concordance Browser

Interactive R/Shiny browser for the Molecular Concordance Index cardiomyopathy project.

## Main input table

`data/mci_concordance_table.csv`

This table is derived from:

`results/mci_scores/TASK6_FINAL_HCM_MCI_WITH_GTEx_AND_BOOTSTRAP_CI.csv`

## Per-cohort DE table

`data/mci_per_cohort_de_table.csv`

This table is derived from:

`results/mci_scores/TASK4_CODE23_HCM_only_clean_DE_input_used_for_MCI.csv`

## Features

1. Gene-level MCI query
2. Per-cohort visual concordance profile
3. GTEx baseline overlay
4. Stratum browser
5. Full concordance table
6. CSV download

## Local run

In R:

setwd('shiny_app/mci_concordance_browser')
shiny::runApp()

## Deployment

This app is intended for shinyapps.io deployment after local testing.

## Interpretation note

MCI measures transcript-level reproducibility across disease cohorts.
GTEx adjustment adds a normal left-ventricle baseline comparison.
GTEx low-confidence status does not mean a gene is biologically irrelevant;
it means disease-associated variability did not exceed normal GTEx left-ventricle variability in this benchmark.