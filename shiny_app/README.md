# MCI Cardiomyopathy Resource Browser

This Shiny app provides an interactive browser for the MCI cardiomyopathy resource.

## Input table

../results/resource_tables/MCI_BROWSER_READY_DISPLAY_TABLE_v0_2.csv

## Main features

- Gene-level search
- MCI, adjusted MCI, tier, and bootstrap interval display
- GTEx low-confidence flag display
- DCM generalization fields where available
- Full searchable and downloadable resource table
- Tier distribution plot
- MCI versus GTEx variability-ratio plot

## Local run

From the repository root, run in R:

shiny::runApp('shiny_app')

Required R packages:

install.packages(c('shiny', 'DT', 'ggplot2'))

## Interpretation boundary

This app is an evidence-auditing browser. It does not prove clinical actionability, therapeutic validity, or disease causality.
