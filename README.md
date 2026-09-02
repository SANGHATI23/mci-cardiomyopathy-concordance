# Molecular Concordance Index (MCI) for Cardiomyopathy

A reproducible computational framework for evaluating the **cross-cohort transcriptomic concordance of ClinVar-annotated cardiomyopathy genes** across public hypertrophic cardiomyopathy (HCM) and dilated cardiomyopathy (DCM) datasets.

The project asks a simple but important question:

> **When a gene is clinically associated with cardiomyopathy, how reproducible is its disease-associated transcriptomic signal across independent human cohorts?**

Rather than treating clinical pathogenicity annotation and transcriptomic reproducibility as the same form of evidence, the **Molecular Concordance Index (MCI)** evaluates whether gene-expression evidence is directionally consistent, similar in magnitude, and statistically reproducible across cohorts.

---

## Why this project?

Genes associated with inherited cardiomyopathy can have strong clinical or genetic evidence while showing heterogeneous transcriptional behavior across human heart datasets.

Differences may arise from:

* cohort composition;
* disease stage;
* tissue sampling;
* treatment exposure;
* technical platform;
* biological heterogeneity;
* sample size;
* baseline expression variability.

Therefore, failure to reproduce a transcriptomic signal does **not** imply that a clinically pathogenic gene is biologically unimportant.

Instead, this project separates several evidence layers that are often conflated:

1. **Clinical pathogenicity annotation**
2. **Cross-cohort transcriptomic reproducibility**
3. **Normal left-ventricle expression variability**
4. **Independent disease-context or external convergence evidence**

The goal is to make those evidence layers transparent, reproducible, and queryable.

---

# Molecular Concordance Index

For gene \(g\):

```text
MCI_g = 0.40 × D_g + 0.35 × S_g + 0.25 × R_g
```

where:

* **D_g — Direction agreement**
  Measures whether independent cohorts show concordant up- or down-regulation.

* **S_g — Effect-size consistency**
  Measures whether estimated disease-associated effects are similar across cohorts.

* **R_g — Statistical reproducibility**
  Measures repeated statistical support across independent cohorts.

The primary weighting scheme therefore gives the largest contribution to direction agreement, followed by effect-size consistency and statistical reproducibility.

### MCI evidence tiers

| MCI score                | Interpretation               |
| ------------------------ | ---------------------------- |
| **≥ 0.70**               | HIGH                         |
| **0.45–0.69**            | MODERATE                     |
| **< 0.45**               | UNSTABLE                     |
| **< 2 eligible cohorts** | INSUFFICIENT_COHORT_COVERAGE |

These tiers describe **transcriptomic evidence stability**, not clinical pathogenicity.

---

# Study workflow

```text
Clinically annotated cardiomyopathy genes
                    ↓
       Eligible human HCM cohorts
                    ↓
     Cohort-specific differential
          expression estimates
                    ↓
       Cross-cohort concordance
                    ↓
    Direction + Effect consistency
      + Statistical reproducibility
                    ↓
          Molecular Concordance
              Index (MCI)
                    ↓
       GTEx normal-LV baseline
             adjustment
                    ↓
      Uncertainty / robustness
             evaluation
                    ↓
       DCM generalization and
          external evidence
                    ↓
     Benchmarking and sensitivity
              analysis
                    ↓
     Queryable cardiomyopathy
           evidence resource
```

---

# Current resource

The current HCM-primary resource contains **49 ClinVar-annotated HCM genes**.

| Evidence category            |  Genes |
| ---------------------------- | -----: |
| HIGH MCI                     | **18** |
| MODERATE MCI                 |  **8** |
| UNSTABLE MCI                 | **21** |
| Insufficient cohort coverage |  **2** |

Additional evidence layers currently include:

* **47 genes** with disease-to-GTEx variability ratios;
* **46 genes** flagged as GTEx-low-confidence under the current rule;
* **3 genes** not flagged as GTEx-low-confidence;
* **22 genes** with attached DCM generalization evidence.

The primary master resource is available at:

```text
results/resource_tables/
MCI_CARDIOMYOPATHY_RESOURCE_MASTER_TABLE_v0_2.csv
```

The corresponding data dictionary is available at:

```text
results/resource_tables/
MCI_CARDIOMYOPATHY_RESOURCE_DATA_DICTIONARY_v0_2.csv
```

For additional interpretation guidance, see [`README_RESOURCE.md`](README_RESOURCE.md).

---

# GTEx baseline adjustment

Cross-cohort disease consistency alone does not establish that the observed variation exceeds normal biological variation.

The GTEx layer therefore compares disease-associated variability with normal human left-ventricle expression variability:

```text
Disease-to-GTEx ratio = σ_disease / σ_GTEx
```

Genes with:

```text
ratio ≤ 1.0
```

are flagged as **GTEx-low-confidence**.

This flag should **not** be interpreted as evidence that the gene is clinically irrelevant.

Instead, it means that, within the current bulk-transcriptomic data, disease-associated variability does not clearly exceed normal left-ventricle baseline variability.

Such genes may require additional evidence from:

* genotype-stratified cohorts;
* single-cell or single-nucleus RNA sequencing;
* proteomics;
* spatial transcriptomics;
* functional assays;
* independent disease cohorts.

---

# Analysis notebooks

The project is organized as a sequential five-notebook workflow.

## 1. Project Setup and Data Acquisition

[`01_MCI_Project_Setup_and_Data_Acquisition.ipynb`](notebooks/01_MCI_Project_Setup_and_Data_Acquisition.ipynb)

Establishes the reproducible foundation of the project.

Main functions include:

* project configuration;
* dataset manifest construction;
* cohort eligibility checks;
* target-gene definition;
* acquisition and organization of source data;
* cohort-specific differential-expression preparation;
* generation of the initial MCI analysis inputs.

This notebook establishes **what data and genes are eligible for analysis and why**.

---

## 2. GTEx Baseline Adjustment

[`02_GTEx_Baseline_Adjustment.ipynb`](notebooks/02_GTEx_Baseline_Adjustment.ipynb)

Adds a normal-tissue reference layer using GTEx left-ventricle expression data.

The analysis asks:

> Is observed disease-associated variability greater than the variability expected in normal human left ventricular tissue?

This helps distinguish apparent disease-associated transcriptomic instability from ordinary baseline expression variability.

---

## 3. DCM Generalization and External Validation

[`03_DCM_Generalization_and_External_Validation.ipynb`](notebooks/03_DCM_Generalization_and_External_Validation.ipynb)

Evaluates whether HCM-derived transcriptomic evidence extends to an additional cardiomyopathy context.

The notebook investigates:

* DCM expression evidence;
* direction of disease-associated effects;
* cross-phenotype agreement;
* disease-context specificity;
* external evidence portability.

This layer helps distinguish genes showing potentially broader cardiomyopathy transcriptomic behavior from genes whose expression evidence appears more phenotype-specific.

---

## 4. MCI Resource and Browser

[`04_MCI_Resource_and_Browser.ipynb`](notebooks/04_MCI_Resource_and_Browser.ipynb)

Transforms the analytical outputs into a structured gene-level evidence resource.

The notebook integrates:

* MCI scores;
* evidence tiers;
* component scores;
* GTEx baseline information;
* uncertainty information;
* DCM generalization evidence;
* provenance and interpretation fields.

Its purpose is to move the project from a collection of analyses to a **queryable evidence-auditing resource**.

Browser-related files are maintained under:

```text
shiny_app/
```

---

## 5. Benchmarking and Sensitivity Analysis

[`05_MCI_Benchmarking_and_Sensitivity_Analysis.ipynb`](notebooks/05_MCI_Benchmarking_and_Sensitivity_Analysis.ipynb)

Evaluates whether the conclusions of MCI depend strongly on its exact formulation.

The primary MCI is compared with several alternative evidence-integration approaches:

1. Equal-weight MCI components
2. Majority direction voting
3. Fisher combined significance
4. Direction-aware Stouffer combined significance
5. DerSimonian-Laird random-effects meta-analysis

Sensitivity analyses also evaluate:

* plausible alternative MCI weights;
* deliberately component-heavy weight schemes;
* leave-one-component-out scoring;
* HIGH/MODERATE threshold changes;
* leave-one-HCM-cohort-out recomputation;
* held-out expression evidence;
* GWAS Catalog convergence.

---

# Benchmarking results

Benchmarking was performed among **47 genes with valid primary MCI scores**.

### Agreement with alternative methods

The equal-weight version of MCI produced nearly identical rankings to the primary weighting scheme:

```text
Spearman ρ = 0.9990
Kendall τ  = 0.9906
```

Only **2 of 47 genes** changed evidence tier, and all **18 original HIGH genes** remained within the equal-weight top-18 set.

Random-effects meta-analysis also showed strong agreement:

```text
Spearman ρ = 0.9427
Kendall τ  = 0.7922
```

and retained **17 of 18** original HIGH genes.

Direction voting showed substantial rank agreement but provided much less discrimination because many genes received identical scores.

Combined-significance approaches such as Fisher and Stouffer produced greater ranking differences because they primarily capture accumulated statistical significance rather than the multidimensional combination of direction, magnitude, and repeated statistical evidence used by MCI.

---

# Weight sensitivity

Three plausible alternative component-weight schemes produced:

```text
Spearman correlations = 0.9973–0.9985
```

relative to the primary MCI.

Only **1–2 of 47 genes** changed tier under these alternative weight schemes.

This indicates that the exact primary coefficients:

```text
0.40 / 0.35 / 0.25
```

are **not the principal driver of the overall gene ranking**.

More substantial changes occurred only when entire MCI components were deliberately removed, supporting the interpretation that MCI behaves as a multidimensional measure rather than as a single-component score.

---

# Threshold sensitivity

Twenty-five combinations of tier thresholds were evaluated:

```text
MODERATE threshold: 0.400–0.500
HIGH threshold:     0.650–0.750
```

Across this grid:

* **39 of 47 genes (83.0%)** retained the same evidence tier under every tested threshold combination.
* Only **8 genes** changed tier anywhere in the grid.
* Changes were concentrated near the original decision boundaries.

The overall resource therefore shows relatively limited sensitivity to reasonable changes in tier thresholds.

---

# Cohort sensitivity

Leave-one-cohort-out analysis revealed a more important source of uncertainty.

Removing individual HCM cohorts produced larger changes than changing MCI weights or tier thresholds.

This indicates that some gene-level classifications remain **cohort dependent**, particularly because removing one cohort from a three-cohort analysis leaves only two independent disease estimates.

Therefore, the project explicitly reports cohort dependence rather than claiming universal cohort-invariant MCI classifications.

This is an important interpretation boundary of the resource.

---

# External evidence

The benchmarking workflow also evaluates MCI against independent evidence sources.

A strict held-out expression comparison produced strong exploratory discrimination:

```text
AUROC            = 0.9884
Average precision = 0.8333
```

However, only **two genes** met the strict held-out positive definition, so this result should be interpreted cautiously.

Directional held-out evidence and GWAS convergence did **not** demonstrate definitive statistical separation.

Accordingly, these analyses are treated as exploratory external evidence rather than proof of predictive validity.

---

# What MCI is — and what it is not

### MCI is intended to be:

* a transparent transcriptomic evidence-auditing framework;
* a measure of cross-cohort gene-expression reproducibility;
* a way to distinguish stable from heterogeneous transcriptomic evidence;
* a framework for integrating disease, normal-tissue, and external evidence layers;
* a resource for research prioritization and translational evidence review.

### MCI is not intended to:

* redefine ClinVar pathogenicity;
* prove that a gene causes cardiomyopathy;
* establish clinical actionability;
* predict patient-level disease risk;
* demonstrate therapeutic efficacy;
* replace genetic or functional evidence;
* claim universal superiority over meta-analysis;
* imply that tier assignments are invariant to cohort composition.

A gene may therefore be **clinically pathogenic yet transcriptomically unstable**, and that distinction is intentional.

---

# Repository structure

```text
mci-cardiomyopathy-concordance/
│
├── data/
│   └── Project input and processed data assets
│
├── metadata/
│   └── Dataset, cohort, and analysis metadata
│
├── notebooks/
│   ├── 01_MCI_Project_Setup_and_Data_Acquisition.ipynb
│   ├── 02_GTEx_Baseline_Adjustment.ipynb
│   ├── 03_DCM_Generalization_and_External_Validation.ipynb
│   ├── 04_MCI_Resource_and_Browser.ipynb
│   └── 05_MCI_Benchmarking_and_Sensitivity_Analysis.ipynb
│
├── scripts/
│   └── Reusable analysis and resource-generation code
│
├── results/
│   ├── differential_expression/
│   ├── figures/
│   ├── logs/
│   ├── mci_scores/
│   ├── quality_control/
│   ├── resource_tables/
│   └── task11_benchmarking/
│
├── shiny_app/
│   └── Interactive-resource/browser components
│
├── manuscript/
│   └── Resource-paper materials
│
├── docs/
│   └── Interpretation and resource-framing documentation
│
├── project_config.json
├── requirements-task11.txt
├── environment-task11.yml
├── README_RESOURCE.md
└── README.md
```

---

# Reproducibility

Clone the repository:

```bash
git clone https://github.com/SANGHATI23/mci-cardiomyopathy-concordance.git
cd mci-cardiomyopathy-concordance
```

For the benchmarking and sensitivity-analysis environment:

```bash
pip install -r requirements-task11.txt
```

or:

```bash
conda env create -f environment-task11.yml
conda activate mci-task11
```

The notebooks are designed to document the analytical workflow from project setup through resource construction and robustness evaluation.

For the benchmarking workflow, run:

```text
notebooks/05_MCI_Benchmarking_and_Sensitivity_Analysis.ipynb
```

from top to bottom.

Benchmarking outputs are written under:

```text
results/task11_benchmarking/
```

Publication-facing tables and figures are available under:

```text
results/task11_benchmarking/publication_ready/
```

The benchmarking input manifest includes file roles, dimensions, sizes, and SHA-256 checksums to support reproducibility auditing.

---

# Resource interpretation

The main value of this project is not assigning a single definitive label to each cardiomyopathy gene.

Instead, the resource allows a user to ask questions such as:

* Is this clinically annotated gene transcriptomically reproducible across HCM cohorts?
* Is its direction of change consistent?
* Are effect sizes stable across studies?
* Is statistical evidence repeatedly observed?
* Does disease-associated variability exceed normal GTEx variability?
* Does the signal extend to DCM?
* Is the MCI classification robust to reasonable parameter choices?
* Does the classification depend strongly on one cohort?
* Is there independent expression or genetic convergence evidence?

The resulting evidence profile is intended to be more informative than a single significance test.

---

# Main interpretation

The analyses support MCI as a **transparent multidimensional framework for auditing transcriptomic evidence across heterogeneous human cardiomyopathy datasets**.

The results suggest that MCI rankings are highly stable to reasonable changes in component weights and relatively stable to changes in tier thresholds.

At the same time, leave-one-cohort-out analyses demonstrate meaningful cohort dependence for some genes.

The appropriate conclusion is therefore not that MCI provides a universal biological truth, but that it provides a reproducible way to make **agreement, disagreement, uncertainty, baseline variability, and evidence portability visible at the gene level**.

---

# Project status

**Current resource version:** `v0.2-resource-reframe-corrected-source`

Current development focuses on:

* cardiomyopathy evidence-resource refinement;
* manuscript/resource-paper development;
* browser/resource presentation;
* additional orthogonal evidence integration;
* future single-cell contextualization.

---

# Citation

A manuscript describing the MCI cardiomyopathy resource is in development.

Until a manuscript or preprint is available, please cite this repository as:

> Basu S. **Molecular Concordance Index: a reproducible transcriptomic resource for evaluating ClinVar-annotated cardiomyopathy genes across human HCM and DCM cohorts.** GitHub repository.

Repository:

https://github.com/SANGHATI23/mci-cardiomyopathy-concordance

---

# Author

**Sanghati Basu**

Health Informatics | Clinical Data Science | Biomedical Informatics

---

# License

This repository is distributed under the terms of the [MIT License](LICENSE).

