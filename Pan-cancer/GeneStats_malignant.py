#!/usr/bin/env python
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse

def gene_expr_stats_malignant(adata_sub):
    mal_mask = adata_sub.obs["celltype"] == "Malignant"
    adata_mal = adata_sub[mal_mask]
    n_cells = adata_mal.n_obs
    X = adata_mal.X
    
    if issparse(X):
        n_expressed = np.asarray((X > 0).sum(axis=0)).flatten()
        
        gene_sum = np.asarray(X.sum(axis=0)).flatten()
        mean_expr = gene_sum / n_cells
        
        X_sq = X.copy()
        X_sq.data **= 2
        mean_sq = np.asarray(X_sq.sum(axis=0)).flatten() / n_cells
        std_expr = np.sqrt(np.maximum(mean_sq - mean_expr ** 2, 0))
    else:
        n_expressed = (X > 0).sum(axis=0)
        mean_expr = X.mean(axis=0)
        std_expr = X.std(axis=0)
    
    pct_expressed = n_expressed / n_cells * 100
    
    stats = pd.DataFrame({
        "gene": adata_mal.var_names,
        "n_cells_expressed": n_expressed.astype(int),
        "pct_cells_expressed": pct_expressed,
        "mean_expression": mean_expr,
        "std_expression": std_expr,
    }).set_index("gene")
    
    print(f"  → {stats.shape[0]:,} genes  |  {n_cells:,} malignant cells used")
    return stats

adata = sc.read_h5ad(
    "/mnt/tank/scratch/achesnokova/flamingo/korean_tisch_3ca_data/merged_all_464.h5ad"
)

stats = gene_expr_stats_malignant(adata)

stats.to_csv("/mnt/tank/scratch/achesnokova/flamingo/korean_tisch_3ca_data/malignant_gene_expression_stats.csv")
