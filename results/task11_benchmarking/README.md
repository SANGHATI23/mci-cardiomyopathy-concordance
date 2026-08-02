# Task 11: MCI Benchmarking and Sensitivity Analysis

Generated: 2026-08-02

## Purpose

Task 11 evaluates whether the Molecular Concordance Index (MCI) behaves robustly relative to simpler benchmark methods and whether rankings or tiers depend materially on component weights, decision thresholds, or individual HCM cohorts.

## Analysis notebook

`notebooks/TASK11_MCI_BENCHMARKING_AND_SENSITIVITY_ANALYSIS.ipynb`

## Benchmark methods

The locked primary MCI was compared with:

1. Equal-weight composite of D_g, S_g, and R_g.
2. Majority direction-vote score.
3. Fisher combined significance.
4. Direction-aware Stouffer combined significance.
5. DerSimonian-Laird random-effects meta-analysis.

## Sensitivity analyses

The notebook evaluates:

- equal and plausible alternative component weights;
- deliberately component-heavy adversarial weights;
- leave-one-component-out scores with renormalization;
- MODERATE thresholds from 0.400 to 0.500;
- HIGH thresholds from 0.650 to 0.750;
- leave-one-HCM-cohort-out recomputation;
- held-out GSE160997 expression evidence;
- GWAS Catalog HCM convergence.

## Exact inputs

The complete input manifest, including file roles, dimensions, sizes, and SHA-256 checksums, is located at:

`results/task11_benchmarking/TASK11_INPUT_DATA_MANIFEST.csv`

## Environment setup

Using pip:

```bash
pip install -r requirements-task11.txt
```

Using Conda:

```bash
conda env create -f environment-task11.yml
conda activate mci-task11
```

## Run instructions

```bash
git clone https://github.com/SANGHATI23/mci-cardiomyopathy-concordance.git
cd mci-cardiomyopathy-concordance
pip install -r requirements-task11.txt
jupyter notebook notebooks/TASK11_MCI_BENCHMARKING_AND_SENSITIVITY_ANALYSIS.ipynb
```

Run all notebook cells sequentially from top to bottom. All generated files are written beneath `results/task11_benchmarking/`.

## Deterministic execution

Task 11 contains no new stochastic resampling. Benchmark calculations, rank correlations, weight tests, threshold-grid tests, leave-one-cohort-out analyses, external-evidence comparisons, and figures are deterministic.

The earlier bootstrap results are consumed as locked inputs and are not regenerated in this notebook.

Additional conventions:

- descending ranks use average ranks for tied values;
- top-N comparisons include every gene tied at the cutoff;
- p-values are bounded at 1e-300 before -log10 transformation;
- random-effects pooling uses the DerSimonian-Laird between-study variance estimator.

## Expected outputs

The notebook generates:

- gene-level scores for all benchmark methods;
- rank-agreement and displacement tables;
- equal-weight tier comparisons;
- weight-sensitivity summaries and gene-level audits;
- threshold-grid summaries and full tier matrices;
- leave-one-cohort-out summaries and gene-level audits;
- held-out and GWAS benchmarking results;
- manuscript-facing compact and full tables;
- five figures exported as PNG, PDF, and SVG;
- input and publication-output manifests.

Publication-ready files are located at:

`results/task11_benchmarking/publication_ready/`

## Interpretation boundary

Task 11 evaluates reproducibility, ranking agreement, and evidence-layer behavior. It does not establish MCI as a clinically validated predictor and does not demonstrate statistical superiority over every simpler method.

The intended use remains transparent transcriptomic evidence auditing and translational triage.
