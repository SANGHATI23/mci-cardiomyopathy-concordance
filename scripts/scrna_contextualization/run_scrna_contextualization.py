
"""
scRNA-seq contextualization module for the MCI cardiomyopathy resource.

This script is a scaffold. It is designed for public human cardiac
single-cell or single-nucleus RNA-seq datasets loaded as AnnData `.h5ad`
files.

Expected input:
    1. AnnData file with genes in adata.var_names
    2. Cell metadata column containing cell-type labels
    3. MCI resource master table CSV

Expected output:
    Gene-level cell-type contextualization table.

Example:
    python scripts/scrna_contextualization/run_scrna_contextualization.py \
        --h5ad data/external/scrna/example_human_heart.h5ad \
        --mci_table results/resource_tables/MCI_CARDIOMYOPATHY_RESOURCE_MASTER_TABLE_v0_2.csv \
        --celltype_col cell_type \
        --output results/resource_tables/MCI_scRNA_celltype_context_v0_1.csv
"""

import argparse
from pathlib import Path
import pandas as pd
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, help="Path to AnnData .h5ad file.")
    parser.add_argument("--mci_table", required=True, help="Path to MCI resource master table.")
    parser.add_argument("--celltype_col", required=True, help="Cell-type column in adata.obs.")
    parser.add_argument("--output", required=True, help="Output CSV path.")
    return parser.parse_args()


def main():
    args = parse_args()

    try:
        import scanpy as sc
    except ImportError as exc:
        raise ImportError(
            "scanpy is required for this module. Install with: pip install scanpy anndata"
        ) from exc

    h5ad_path = Path(args.h5ad)
    mci_path = Path(args.mci_table)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(h5ad_path)
    mci = pd.read_csv(mci_path)

    if args.celltype_col not in adata.obs.columns:
        raise ValueError(f"Cell-type column not found in adata.obs: {args.celltype_col}")

    if "gene_symbol" not in mci.columns:
        raise ValueError("MCI table must contain gene_symbol column.")

    # Standardize gene symbols
    mci["gene_symbol"] = mci["gene_symbol"].astype(str).str.upper().str.strip()
    adata.var_names = adata.var_names.astype(str).str.upper().str.strip()

    target_genes = sorted(set(mci["gene_symbol"]) & set(adata.var_names))
    missing_genes = sorted(set(mci["gene_symbol"]) - set(adata.var_names))

    if len(target_genes) == 0:
        raise ValueError("No MCI target genes found in AnnData var_names.")

    # Subset target genes
    adata_sub = adata[:, target_genes].copy()

    # Convert to dense only for small target-gene subset
    X = adata_sub.X
    if hasattr(X, "toarray"):
        X = X.toarray()

    expr = pd.DataFrame(X, columns=target_genes, index=adata_sub.obs_names)
    expr[args.celltype_col] = adata_sub.obs[args.celltype_col].astype(str).values

    # Mean expression by cell type
    means = expr.groupby(args.celltype_col)[target_genes].mean().T
    means.index.name = "gene_symbol"
    means = means.reset_index()

    # Long format for specificity computation
    long_df = means.melt(id_vars="gene_symbol", var_name="cell_type", value_name="mean_expression")

    # Dominant cell type and specificity
    dominant = (
        long_df.sort_values(["gene_symbol", "mean_expression"], ascending=[True, False])
        .groupby("gene_symbol")
        .head(1)
        .rename(columns={"cell_type": "dominant_cell_type", "mean_expression": "dominant_celltype_mean_expression"})
    )

    total_expr = long_df.groupby("gene_symbol")["mean_expression"].sum().reset_index()
    total_expr = total_expr.rename(columns={"mean_expression": "sum_mean_expression_across_celltypes"})

    dominant = dominant.merge(total_expr, on="gene_symbol", how="left")
    dominant["cell_type_specificity_score"] = (
        dominant["dominant_celltype_mean_expression"] /
        dominant["sum_mean_expression_across_celltypes"].replace(0, np.nan)
    )

    # Wide cell-type table with safe column names
    means_wide = means.copy()
    safe_cols = []
    for col in means_wide.columns:
        if col == "gene_symbol":
            safe_cols.append(col)
        else:
            safe = str(col).lower().replace(" ", "_").replace("/", "_").replace("-", "_")
            safe_cols.append(f"mean_expression_{safe}")
    means_wide.columns = safe_cols

    out = mci.merge(means_wide, on="gene_symbol", how="left")
    out = out.merge(
        dominant[[
            "gene_symbol",
            "dominant_cell_type",
            "dominant_celltype_mean_expression",
            "cell_type_specificity_score"
        ]],
        on="gene_symbol",
        how="left"
    )

    out["matched_in_scrna_dataset"] = out["gene_symbol"].isin(target_genes)
    out["scrna_missing_reason"] = np.where(
        out["matched_in_scrna_dataset"],
        "",
        "gene_symbol_not_found_in_scrna_var_names"
    )

    out.to_csv(output_path, index=False)

    summary = {
        "h5ad": str(h5ad_path),
        "mci_table": str(mci_path),
        "output": str(output_path),
        "celltype_col": args.celltype_col,
        "n_mci_genes": int(mci["gene_symbol"].nunique()),
        "n_matched_genes": int(len(target_genes)),
        "n_missing_genes": int(len(missing_genes)),
        "missing_genes": missing_genes[:100],
        "n_cells": int(adata.n_obs),
        "n_genes_in_scrna": int(adata.n_vars),
        "cell_type_counts": adata.obs[args.celltype_col].astype(str).value_counts().to_dict()
    }

    summary_path = output_path.with_suffix(".summary.json")
    summary_path.write_text(pd.Series(summary).to_json(indent=2))

    print("scRNA contextualization complete.")
    print("Output:", output_path)
    print("Summary:", summary_path)


if __name__ == "__main__":
    main()
