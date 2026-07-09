
# Molecular Concordance Index: a reproducible transcriptomic resource for evaluating ClinVar-annotated cardiomyopathy genes across human HCM and DCM cohorts

Sanghati Basu and collaborators

## Abstract

Clinical pathogenicity annotation is essential for cardiomyopathy genetics, but it does not establish whether a gene shows reproducible transcript-level perturbation across independent human disease cohorts. This creates a translational informatics gap: genes may be clinically important while remaining unstable, cohort-specific, or indistinguishable from normal-tissue variability at the bulk-transcriptome level. We present the Molecular Concordance Index, or MCI, a reproducible resource and scoring framework for evaluating transcriptomic concordance of ClinVar-annotated hypertrophic cardiomyopathy and dilated cardiomyopathy genes across public human heart datasets.

MCI integrates three evidence layers: direction agreement across cohorts, effect-size consistency, and statistical reproducibility. We applied the framework to public HCM cohorts, extended the analysis to DCM generalization datasets, added 1,000-iteration cohort-level bootstrap confidence intervals, and benchmarked disease-associated variability against GTEx v8 Heart - Left Ventricle expression variability. The resulting resource assigns each gene a concordance score, bootstrap interval, tier classification, GTEx low-confidence flag, and per-cohort differential-expression evidence.

Across 49 HCM genes, the resource identified heterogeneous transcript-level behavior: 18 genes were classified as high concordance, 8 as moderate, 21 as unstable, and 2 as insufficient coverage. GTEx benchmarking showed that 46 of 49 genes were flagged as GTEx-low-confidence, indicating that disease-associated variability often did not exceed normal left-ventricle baseline variability in the current bulk data. DCM generalization entries were available for 22 genes. Pre-specified mechanism-stratified, DCM generalization, held-out replication, and GWAS convergence analyses were directionally informative but did not provide definitive statistically significant discovery claims.

The MCI resource therefore provides a transparent evidence-auditing layer for cardiomyopathy target evaluation rather than a binary discovery test. It enables researchers to distinguish clinically annotated genes with reproducible transcriptomic support from genes whose disease signal is unstable, baseline-confounded, disease-context-specific, or underpowered. All resource tables, figures, and browser-ready outputs are designed for public reuse and extension, including future single-cell and genotype-stratified modules.

## Keywords

cardiomyopathy; ClinVar; transcriptomics; target validation; molecular concordance; GTEx; hypertrophic cardiomyopathy; dilated cardiomyopathy; database resource; biomedical informatics

## 1. Introduction

Inherited cardiomyopathies are genetically heterogeneous disorders in which clinical interpretation often relies on curated pathogenic and likely pathogenic variant annotations. ClinVar and related clinical-genomic resources are therefore essential for identifying genes implicated in hypertrophic cardiomyopathy and dilated cardiomyopathy. However, clinical pathogenicity annotation and transcriptomic reproducibility are not the same evidence layer. A gene may be clinically causal through protein structure, sarcomere mechanics, splicing, ion handling, or cell-type-specific mechanisms without showing reproducible bulk-transcript abundance change across independent disease cohorts.

This distinction matters for translational informatics. Public human heart transcriptomic datasets are increasingly used to prioritize biomarkers, nominate targets, and justify downstream validation. Yet there is no standard resource that asks whether ClinVar-annotated cardiomyopathy genes show reproducible direction, effect size, and statistical evidence across independent human heart cohorts. There is also limited separation between disease-cohort reproducibility and baseline variability in normal human left ventricle.

The Molecular Concordance Index resource addresses this gap. Instead of treating one differential-expression result as sufficient evidence, MCI evaluates whether a gene behaves consistently across cohorts. It combines direction agreement, effect-size consistency, and statistical reproducibility into a gene-level concordance score. The resource further adds bootstrap uncertainty, GTEx v8 Heart - Left Ventricle baseline benchmarking, DCM generalization, held-out validation, GWAS convergence, and browser-ready outputs.

The purpose of this manuscript is therefore not to claim a single definitive cardiomyopathy mechanism discovery. The purpose is to present a reusable evidence-auditing resource. The core claim is that clinical pathogenicity annotation, bulk transcriptomic concordance, normal-tissue expression variability, and independent genetic association are related but non-equivalent layers of evidence.

## 2. Results

### 2.1 Resource construction and cohort coverage

The MCI resource was constructed from public human heart transcriptomic cohorts, ClinVar-annotated cardiomyopathy gene sets, and GTEx v8 Heart - Left Ventricle baseline expression. The current HCM-primary resource contains 49 ClinVar-annotated genes with MCI scoring or coverage status.

The primary resource table is:

`results/resource_tables/MCI_CARDIOMYOPATHY_RESOURCE_MASTER_TABLE_v0_2.csv`

A data dictionary is provided at:

`results/resource_tables/MCI_CARDIOMYOPATHY_RESOURCE_DATA_DICTIONARY_v0_2.csv`

Each row includes gene symbol, disease group, pre-specified stratum where available, MCI, adjusted MCI, bootstrap interval, tier assignment, GTEx baseline fields, DCM generalization fields where available, and a plain-language interpretation field for browser display.

### 2.2 Primary HCM MCI resource

Across 49 HCM genes, 18 were classified as HIGH, 8 as MODERATE, 21 as UNSTABLE, and 2 as INSUFFICIENT_COHORT_COVERAGE. The median MCI was 0.584, and the mean MCI was 0.537.

This distribution supports the main resource premise: clinically annotated cardiomyopathy genes do not behave uniformly at the bulk-transcriptomic level. Some genes show strong cross-cohort concordance, whereas others remain unstable or insufficiently covered.

#### Top high-concordance HCM genes in the current resource

| gene_symbol   |      MCI |   Adj_MCI | MCI_tier   | GTEx_low_confidence_flag   |
|:--------------|---------:|----------:|:-----------|:---------------------------|
| MYH6          | 0.96806  |  1.18436  | HIGH       | True                       |
| TNNI3         | 0.961721 |  1.05051  | HIGH       | True                       |
| GLA           | 0.960029 |  1.18518  | HIGH       | True                       |
| ACTN2         | 0.934792 |  1.00024  | HIGH       | True                       |
| TMEM43        | 0.929496 |  1.21375  | HIGH       | True                       |
| TINF2         | 0.904054 |  1.19951  | HIGH       | True                       |
| KCNH2         | 0.881352 |  1.12804  | HIGH       | True                       |
| BAG3          | 0.828181 |  0.976573 | HIGH       | True                       |
| TRIM63        | 0.819793 |  0.832012 | HIGH       | True                       |
| RBM20         | 0.803574 |  1.07943  | HIGH       | True                       |

#### Most unstable HCM genes in the current resource

| gene_symbol   |      MCI | MCI_tier   | GTEx_low_confidence_flag   |
|:--------------|---------:|:-----------|:---------------------------|
| CEP85L        | 0        | UNSTABLE   | True                       |
| DES           | 0.133333 | UNSTABLE   | True                       |
| JPH2          | 0.133333 | UNSTABLE   | True                       |
| CACNA1C       | 0.133333 | UNSTABLE   | True                       |
| MYBPC3        | 0.133333 | UNSTABLE   | True                       |
| MYO6          | 0.133333 | UNSTABLE   | True                       |
| TTN           | 0.133333 | UNSTABLE   | True                       |
| TNNT2         | 0.214271 | UNSTABLE   | True                       |
| MTO1          | 0.216667 | UNSTABLE   | True                       |
| NEXN          | 0.216667 | UNSTABLE   | True                       |

### 2.3 Bootstrap uncertainty and tier stability

The resource includes cohort-level bootstrap uncertainty estimates. Bootstrap intervals are intentionally reported because the current public HCM evidence base contains a limited number of eligible cohorts. The intervals should therefore be interpreted as transparency measures rather than as high-powered uncertainty estimates.

The bootstrap layer supports a database-resource interpretation: users can inspect not only the point MCI tier, but also whether the tier is stable under cohort-level resampling and whether the confidence interval crosses tier thresholds.

### 2.4 GTEx baseline benchmarking

GTEx v8 Heart - Left Ventricle expression was used as a normal baseline benchmark. For each gene, disease-associated variability was compared against normal left-ventricle variability using:

`sigma_disease_to_GTEx_ratio = sigma_disease / sigma_GTEx`

In the current resource, 47 genes had available disease-to-GTEx variability ratios. 46 genes were flagged as GTEx-low-confidence, while 3 were not flagged as low-confidence.

A GTEx-low-confidence flag does not mean that a gene is biologically irrelevant. It means that, in the current bulk RNA-seq and microarray evidence, disease-associated variability does not exceed normal GTEx left-ventricle baseline variability. This is especially important for target interpretation because a gene can be clinically pathogenic and transcript-concordant while still requiring orthogonal support from single-cell, proteomic, genotype-stratified, or functional evidence.

### 2.5 Mechanism-stratified use case

A pre-specified mechanism-stratified analysis compared sarcomeric and non-sarcomeric HCM genes. This analysis should be interpreted as a demonstration of how biological strata can be evaluated within the MCI resource, not as the central discovery claim of the paper. The observed pattern was directionally informative but not statistically definitive.

This result is still useful in a resource-paper framework because it identifies where additional cohorts, genotype-stratified expression, proteomics, or single-cell analysis are needed before making strong mechanism-level claims.

### 2.6 DCM generalization module

The DCM layer evaluates whether concordance patterns observed in HCM extend into DCM datasets. In the current master table, 22 genes have attached DCM generalization entries.

The DCM module is not framed as proof that HCM and DCM transcriptomic concordance behave identically. Instead, it is a disease-context test. Where DCM patterns do not mirror HCM patterns, the resource highlights disease specificity rather than suppressing discordant evidence.

### 2.7 Held-out and GWAS validation modules

The resource includes held-out and GWAS convergence modules as independent evidence layers. These modules should be reported transparently even when they are directionally supportive but not statistically significant. Their value is to show whether MCI tiers align with unseen cohort behavior or inherited genetic association evidence.

In the resource-paper framing, these analyses are not required to produce a binary significant result. They define the current boundary of evidence and clarify where transcriptomic concordance, replication behavior, and GWAS signals agree or diverge.

### 2.8 Browser-ready output

The master table is designed for browser deployment. A gene-level browser should allow users to search a gene and retrieve:

- MCI and adjusted MCI
- MCI tier and bootstrap-majority tier
- bootstrap confidence interval
- GTEx low-confidence flag
- disease-to-GTEx variability ratio
- DCM generalization entry where available
- plain-language interpretation
- downloadable full table

This makes the project suitable for database/resource manuscript tracks.

### 2.9 Planned single-cell contextualization module

A scaffolded scRNA-seq contextualization module is included in:

`scripts/scrna_contextualization/`

This module is designed to annotate MCI genes by cardiac cell-type expression using public human single-cell or single-nucleus RNA-seq datasets. The goal is to determine whether unstable bulk MCI behavior may reflect cell-type dilution or cell-type specificity rather than absence of biological relevance.

## 3. Methods

### 3.1 Dataset acquisition

Public human heart transcriptomic datasets were used for HCM and DCM concordance scoring. The resource uses existing public data and does not require new sample collection.

### 3.2 ClinVar cardiomyopathy gene universe

ClinVar pathogenic and likely pathogenic cardiomyopathy annotations were used to define the gene universe. Genes were grouped into disease and mechanism strata where possible.

### 3.3 Differential-expression processing

Each cohort was analyzed independently to estimate disease-versus-control differential expression. The output for each cohort included gene symbol, log2 fold change, standard error where available, nominal p-value, and FDR-adjusted p-value.

### 3.4 Harmonization and batch-aware processing

Gene identifiers were standardized to HGNC-style gene symbols. Cohort-level outputs were harmonized into a common schema before MCI scoring. Batch-aware processing and harmonization steps were used where appropriate for expression preprocessing and cross-cohort comparability.

### 3.5 Molecular Concordance Index formula

For gene g, MCI was computed as:

`MCI_g = 0.40 * D_g + 0.35 * S_g + 0.25 * R_g`

where `D_g` is direction agreement, `S_g` is effect-size consistency, and `R_g` is statistical reproducibility.

### 3.6 Tier assignment

Genes were assigned to tiers using pre-specified thresholds:

- HIGH: MCI >= 0.70
- MODERATE: 0.45 <= MCI < 0.70
- UNSTABLE: MCI < 0.45
- INSUFFICIENT_COHORT_COVERAGE: fewer than two eligible cohorts

### 3.7 Bootstrap confidence intervals

Cohort-level bootstrap resampling was performed to estimate MCI uncertainty and tier stability. Bootstrap outputs include 95 percent confidence intervals, tier probabilities, and bootstrap-majority tier labels.

### 3.8 GTEx baseline adjustment

GTEx v8 Heart - Left Ventricle expression was used to calculate normal baseline variability. Disease-associated variability was compared against GTEx variability using `sigma_disease_to_GTEx_ratio`. Genes with ratio less than or equal to 1.0 were flagged as GTEx-low-confidence.

### 3.9 DCM generalization

DCM datasets were processed using the same concordance logic where available. DCM entries were attached to the master resource table as a disease-context generalization layer.

### 3.10 Held-out and GWAS validation

Held-out validation and GWAS convergence were implemented as independent validation modules. These results are interpreted as evidence layers rather than as required binary significance tests.

### 3.11 scRNA-seq contextualization scaffold

A planned scRNA-seq module was scaffolded to support future annotation of MCI genes by cardiac cell type. The module accepts public AnnData files and reports gene-level cell-type expression, dominant cell type, and cell-type specificity metrics.

## 4. Discussion

The MCI cardiomyopathy resource addresses a practical translational informatics problem: clinical pathogenicity annotation does not automatically establish reproducible bulk-transcriptomic evidence. By making transcriptomic concordance, GTEx baseline variability, bootstrap uncertainty, and validation modules visible gene by gene, the resource helps researchers avoid overinterpreting single-cohort differential-expression results.

The main finding is not that one mechanism hypothesis is definitively confirmed. The main finding is that cardiomyopathy genes differ substantially in transcript-level reproducibility and baseline-context confidence. This is exactly why a resource is needed.

The GTEx layer is particularly important. Many high-MCI genes are transcript-concordant but GTEx-low-confidence, meaning their disease-associated variability does not exceed normal left-ventricle baseline variability in the current bulk data. This does not invalidate the genes. It changes the interpretation from unqualified transcript-level target to candidate requiring orthogonal validation.

The non-significant mechanism, DCM, held-out, and GWAS analyses should be interpreted as honest evidence boundaries. In a discovery manuscript, these results would weaken the central claim. In a resource manuscript, they strengthen the credibility of the resource because they show that the framework reports discordance and uncertainty rather than hiding it.

## 5. Limitations

The current resource is limited by the number and size of public cohorts, bulk-tissue composition effects, incomplete genotype stratification, platform differences, and lack of direct wet-lab validation. GTEx provides a useful normal baseline but is not a matched disease-control cohort. DCM generalization remains limited by available datasets and gene coverage. The scRNA-seq module has been scaffolded but still requires selection and execution on a suitable public human cardiac single-cell dataset.

## 6. Conclusion

The MCI cardiomyopathy resource provides a reusable and transparent evidence-auditing layer for ClinVar-annotated cardiomyopathy genes. It distinguishes clinical annotation from transcriptomic reproducibility, normal-tissue variability, disease-context generalization, and independent validation evidence. This resource-paper framing better matches the current evidence than a high-impact discovery claim and provides a stronger foundation for future single-cell, proteomic, genotype-stratified, and large-cohort extensions.

## Data and code availability

Code and resource outputs are available in this repository. The current master resource table is located at:

`results/resource_tables/MCI_CARDIOMYOPATHY_RESOURCE_MASTER_TABLE_v0_2.csv`

The data dictionary is located at:

`results/resource_tables/MCI_CARDIOMYOPATHY_RESOURCE_DATA_DICTIONARY_v0_2.csv`

The scRNA-seq contextualization scaffold is located at:

`scripts/scrna_contextualization/`
