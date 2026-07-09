
# scRNA-seq Contextualization Module

## Purpose

This module is designed to extend the MCI cardiomyopathy resource from bulk-tissue transcriptomic concordance into cell-type contextualization using public human cardiac single-cell or single-nucleus RNA-seq datasets.

The goal is not to replace the bulk MCI score. The goal is to explain whether bulk-transcript instability may reflect cell-type dilution, cell-type specificity, or disease-cell-composition effects.

## Scientific motivation

A gene can be clinically important but appear unstable in bulk RNA-seq for several reasons:

1. The gene may be expressed mainly in cardiomyocytes but diluted by fibroblast, endothelial, smooth muscle, or immune-cell composition.
2. The disease signal may exist only in a specific cell type.
3. The gene may act through protein-level or splicing-level mechanisms rather than bulk abundance.
4. Public bulk cohorts may differ in cellular composition, disease stage, tissue region, or genotype mix.

A single-cell contextualization layer can therefore strengthen the resource-paper framing by showing how MCI genes behave across cardiac cell types.

## Planned input

A public human HCM or DCM single-cell/single-nucleus RNA-seq dataset with:

- expression matrix or AnnData object,
- cell barcode metadata,
- cell-type labels,
- disease status if available,
- human gene symbols.

## Planned output

The planned output is a gene-level annotation table:

- gene_symbol
- cardiomyocyte_mean_expression
- fibroblast_mean_expression
- endothelial_mean_expression
- immune_mean_expression
- smooth_muscle_mean_expression
- dominant_cell_type
- cell_type_specificity_score
- detected_in_cardiomyocytes
- detected_in_fibroblasts
- MCI_tier
- GTEx_low_confidence_flag
- interpretation

## Manuscript role

This module should be presented as a cell-type contextualization layer for the database/resource paper. It should not be overclaimed as wet-lab validation.

Recommended sentence:

To support cell-type interpretation of bulk concordance patterns, we added a planned scRNA-seq contextualization module that annotates MCI genes by cardiac cell-type expression and identifies genes whose apparent bulk instability may reflect cell-type specificity or cellular composition effects.

## Current status

Scaffold created. Dataset selection and execution remain pending.
