
"""
MCI scoring functions for the Molecular Concordance Index project.

Formula:
MCI_g = 0.40 * D_g + 0.35 * S_g + 0.25 * R_g

Where:
D_g = direction agreement
S_g = effect size consistency
R_g = statistical reproducibility
"""

import numpy as np
import pandas as pd


def direction_agreement(log2fc_values):
    """
    D_g = proportion of cohort pairs with concordant log2FC sign.
    """
    signs = np.sign(log2fc_values)

    pair_results = []
    for i in range(len(signs)):
        for j in range(i + 1, len(signs)):
            pair_results.append(1 if signs[i] == signs[j] else 0)

    if len(pair_results) == 0:
        return np.nan

    return float(np.mean(pair_results))


def effect_size_consistency(log2fc_values):
    """
    S_g = 1 - min(CV(log2FC) / 2.0, 1.0)

    Uses absolute mean to avoid sign-related instability.
    """
    values = np.array(log2fc_values, dtype=float)

    mean_abs = abs(np.mean(values))
    sd = np.std(values, ddof=1)

    if mean_abs == 0:
        return 0.0

    cv = sd / mean_abs
    return float(1 - min(cv / 2.0, 1.0))


def statistical_reproducibility(fdr_values, threshold=0.05):
    """
    R_g = proportion of cohorts with FDR < threshold.
    """
    return float(np.mean(np.array(fdr_values) < threshold))


def assign_mci_tier(mci):
    """
    Pre-specified MCI tier thresholds.
    """
    if mci >= 0.70:
        return "High"
    elif mci >= 0.45:
        return "Moderate"
    else:
        return "Unstable"


def compute_mci_for_gene(group):
    """
    Compute MCI components and composite score for one gene.
    Input group must contain: log2FC, FDR
    """
    log2fc_values = group["log2FC"].values
    fdr_values = group["FDR"].values

    D_g = direction_agreement(log2fc_values)
    S_g = effect_size_consistency(log2fc_values)
    R_g = statistical_reproducibility(fdr_values)

    MCI_g = (0.40 * D_g) + (0.35 * S_g) + (0.25 * R_g)

    return pd.Series({
        "n_cohorts": len(group),
        "D_g_direction_agreement": D_g,
        "S_g_effect_size_consistency": S_g,
        "R_g_statistical_reproducibility": R_g,
        "MCI": MCI_g,
        "MCI_tier": assign_mci_tier(MCI_g)
    })


def compute_mci_table(de_results):
    """
    Compute MCI scores for all genes in a differential-expression table.

    Required columns:
    - gene_symbol
    - cohort
    - log2FC
    - FDR
    """
    required_columns = {"gene_symbol", "cohort", "log2FC", "FDR"}

    missing_columns = required_columns - set(de_results.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    mci_results = (
        de_results
        .groupby("gene_symbol", group_keys=False)
        .apply(compute_mci_for_gene)
        .reset_index()
    )

    return mci_results
