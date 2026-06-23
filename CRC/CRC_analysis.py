# srun --time=05:00:00 --mem=400G --cpus-per-task=4 --qos=achesnokova --pty bash
# conda activate atlas

import scanpy as sc
import pandas as pd
import numpy as np
import gc
from scipy.sparse import issparse
import re

def group_cell_types(cell_type_list):
    mapping_rules = [
        # Plasma cells
        (r'Ig[AGM] plasma cell', 'plasma cell'),
        (r'plasmablast', 'plasma cell'),
        
        # B cells
        (r'germinal center B cell', 'B cell'),
        (r'mature B cell', 'B cell'),
        (r'memory B cell', 'B cell'),
        (r'naive B cell', 'B cell'),
        
        # T cells
        (r'CD4-positive, alpha-beta T cell', 'CD4 T cell'),
        (r'naive thymus-derived CD4-positive, alpha-beta T cell', 'CD4 T cell'),
        (r'CD8-positive, alpha-beta T cell', 'CD8 T cell'),
        (r'naive thymus-derived CD8-positive, alpha-beta T cell', 'CD8 T cell'),
        (r'gamma-delta T cell', 'gamma-delta T cell'),
        (r'mature NK T cell', 'NKT cell'),
        (r'regulatory T cell', 'Treg'),
        
        # Monocytes/Macrophages
        (r'classical monocyte', 'monocyte'),
        (r'non-classical monocyte', 'monocyte'),
        (r'macrophage', 'macrophage'),
        
        # Dendritic cells
        (r'CD1c-positive myeloid dendritic cell', 'dendritic cell'),
        (r'mature conventional dendritic cell', 'dendritic cell'),
        (r'plasmacytoid dendritic cell', 'dendritic cell'),
        
        # Epithelial cells (intestinal)
        (r'BEST4\+ intestinal epithelial cell', 'epithelial cell'),
        (r'colon goblet cell', 'epithelial cell'),
        (r'enterocyte of colon', 'epithelial cell'),
        (r'enteroendocrine cell of colon', 'epithelial cell'),
        (r'tuft cell of colon', 'epithelial cell'),
        (r'intestinal crypt stem cell of colon', 'epithelial cell'),
        (r'transit amplifying cell of colon', 'epithelial cell'),
        
        # Endothelial cells
        (r'endothelial cell of artery', 'endothelial cell'),
        (r'endothelial cell of lymphatic vessel', 'endothelial cell'),
        (r'vein endothelial cell', 'endothelial cell'),
        
        # Other immune cells
        (r'innate lymphoid cell', 'innate lymphoid cell'),
        (r'natural killer cell', 'NK cell'),
        (r'mast cell', 'mast cell'),
        (r'neutrophil', 'neutrophil'),
        (r'granulocyte monocyte progenitor cell', 'progenitor cell'),
        
        # Stromal cells
        (r'fibroblast', 'fibroblast'),
        (r'reticular cell', 'fibroblast'),
        (r'pericyte', 'pericyte'),
        (r'Schwann cell', 'Schwann cell'),
        
        # Malignant
        (r'malignant cell', 'malignant cell'),
    ]
    
    grouped = []
    for cell in cell_type_list:
        matched = False
        for pattern, group_name in mapping_rules:
            if re.match(pattern, cell, re.IGNORECASE):
                grouped.append(group_name)
                matched = True
                break
        if not matched:
            grouped.append(cell)
    
    return grouped

def print_stats(adata, label):
    n_cells   = adata.n_obs
    n_samples = adata.obs["sample_id"].nunique()
    n_mal     = (adata.obs["cell_type"] == "malignant cell").sum()
    print(
        f"[{label}] "
        f"cells={n_cells:,}  |  samples={n_samples:,}  |  malignant cells={n_mal:,}"
    )

adata = sc.read_h5ad("final_crc_atlas-adata.h5ad")
print_stats(adata, "Raw")

mask_platform = adata.obs["platform"].isin(["10x 3p", "10x 5p"])
adata._inplace_subset_obs(mask_platform)
gc.collect()
print_stats(adata, "After platform filter (10x 3p / 10x 5p)")

mask_enrich = adata.obs["enrichment_cell_types"].isin(["mixCD45+CD45-", "naive", "CD45+"])
adata._inplace_subset_obs(mask_enrich)
gc.collect()
print_stats(adata, "After enrichment filter (mixCD45+CD45-, naive, CD45+)")

mask_cancer = ~adata.obs["cancer_type"].isin(["normal", np.nan])
adata._inplace_subset_obs(mask_cancer)
gc.collect()
print_stats(adata, "After cancer_type filter (drop normal/NaN)")

adata.obs["cell_type_compact"] = group_cell_types(adata.obs["cell_type"].tolist())

immune_genes = [
    'ENSG00000262418',  # PTPRC
    'ENSG00000167286', # CD3D
    'ENSG00000168685',  # IL7R
    'ENSG00000172116', # CD8B
    'ENSG00000153563',  # CD8A
    'ENSG00000105374',  # NKG7
    'ENSG00000110848',  # CD69
    'ENSG00000173372',  # C1QA
    'ENSG00000085265',  # FCN1
    'ENSG00000170476',  # MZB1
]
present_genes = [g for g in immune_genes if g in adata.var_names]

malignant_obs_before = adata.obs_names[adata.obs["cell_type"] == "malignant cell"].tolist()
 
malignant_mask   = (adata.obs["cell_type"] == "malignant cell").values
malignant_indices = np.where(malignant_mask)[0]
 
removed_immune_obs = []
 
if len(malignant_indices) > 0 and len(present_genes) > 0:
    expr_matrix = adata[malignant_mask, present_genes].X
 
    if issparse(expr_matrix):
        has_immune_expr = np.asarray((expr_matrix > 0).sum(axis=1)).flatten() > 0
    else:
        has_immune_expr = (expr_matrix > 0).any(axis=1)
 
    remove_indices = malignant_indices[has_immune_expr]
    remove_mask    = np.zeros(adata.n_obs, dtype=bool)
    remove_mask[remove_indices] = True
 
    removed_immune_obs = adata.obs_names[remove_mask].tolist()
 
    print(
        f"\n Immune-gene filter: {has_immune_expr.sum():,} / {len(malignant_indices):,} "
        f"malignant cells express ≥1 immune gene → will be removed"
    )
 
    adata._inplace_subset_obs(~remove_mask)
    gc.collect()
 
print_stats(adata, "After immune-gene filter on malignant cells")

malignant_counts = (
    adata.obs[adata.obs["cell_type"] == "malignant cell"]
    .groupby("sample_id")
    .size()
)

samples_to_keep = malignant_counts[malignant_counts >= 100].index
mask_samples = adata.obs["sample_id"].isin(samples_to_keep)
 
adata._inplace_subset_obs(mask_samples)
gc.collect()
 
print_stats(adata, f"After ≥100-malignant-cells filter")

def pseudobulk_malignant(adata_sub):
    mal_mask = adata_sub.obs["cell_type"] == "malignant cell"
    adata_mal = adata_sub[mal_mask]
 
    samples  = adata_mal.obs["sample_id"].values
    unique_s = np.array(sorted(set(samples)))
 
    X = adata_mal.layers["counts"] if "counts" in adata_mal.layers else adata_mal.X
 
    rows = []
    for sid in unique_s:
        idx = np.where(samples == sid)[0]
        block = X[idx, :]
        if issparse(block):
            s = np.asarray(block.sum(axis=0)).flatten()
        else:
            s = block.sum(axis=0)
        rows.append(s)
 
    pb = pd.DataFrame(
        np.vstack(rows),
        index=unique_s,
        columns=adata_mal.var_names,
    )
    pb.index.name = "sample_id"
    return pb
 
# ── 1.  Pseudobulk WITH immune-gene filter ─────────────
 
pb_filtered = pseudobulk_malignant(adata)
pb_filtered.to_csv("pseudobulk_malignant_immuneEnrichrFiltered.csv") 

# ── 3.  Sample metadata table ─────────────────────────────────────────────────
sample_meta = (
    adata.obs
    .drop_duplicates(subset="sample_id")
    .set_index("sample_id")
    .sort_index()
)
 
# Append per-sample cell counts
sample_meta["n_cells_total"] = adata.obs.groupby("sample_id").size()
sample_meta["n_malignant_cells"] = (
    adata.obs[adata.obs["cell_type"] == "malignant cell"]
    .groupby("sample_id")
    .size()
    .reindex(sample_meta.index, fill_value=0)
)
 
sample_meta.to_csv("sample_metadata_2.csv")

# ── 4.  Cell-type × Sample count matrix ───────────────────────────────────────
ct_sample = (
    adata.obs
    .groupby(["cell_type_compact", "sample_id"])
    .size()
    .unstack(fill_value=0)
)
ct_sample.index.name  = "cell_type_compact"
ct_sample.columns.name = "sample_id"
 
ct_sample.to_csv("celltype_sample_counts_ImmuneEnrichrFiltered.csv")

# ── 5.  Gene expression statistics in malignant cells ─────────────────────────

def gene_expr_stats_malignant(adata_sub):
    mal_mask  = adata_sub.obs["cell_type"] == "malignant cell"
    adata_mal = adata_sub[mal_mask]
    n_cells   = adata_mal.n_obs
    X = adata_mal.layers["counts"] if "counts" in adata_mal.layers else adata_mal.X
    if issparse(X):
        # n_cells_expressed: number of non-zero entries per gene
        n_expressed = np.asarray((X > 0).sum(axis=0)).flatten()
        # mean: sum / n_cells  (sum over all cells, including zeros)
        gene_sum    = np.asarray(X.sum(axis=0)).flatten()
        mean_expr   = gene_sum / n_cells
        # std: E[x^2] - E[x]^2, computed without densifying
        X_sq        = X.copy()
        X_sq.data **= 2
        mean_sq     = np.asarray(X_sq.sum(axis=0)).flatten() / n_cells
        std_expr    = np.sqrt(np.maximum(mean_sq - mean_expr ** 2, 0))
    else:
        n_expressed = (X > 0).sum(axis=0)
        mean_expr   = X.mean(axis=0)
        std_expr    = X.std(axis=0)
    pct_expressed = n_expressed / n_cells * 100
    stats = pd.DataFrame({
        "gene":               adata_mal.var_names,
        "n_cells_expressed":  n_expressed.astype(int),
        "pct_cells_expressed": pct_expressed,
        "mean_expression":    mean_expr,
        "std_expression":     std_expr,
    }).set_index("gene")
 
    print(f"  → {stats.shape[0]:,} genes  |  {n_cells:,} malignant cells used")
    return stats

gene_stats_filtered = gene_expr_stats_malignant(adata)
gene_stats_filtered.to_csv("gene_stats_malignant_immuneFiltered_2.csv")


def fix_anndata_definitively(adata):
    obs_df = adata.obs.copy()
    for col in obs_df.columns:
        if hasattr(obs_df[col].dtype, 'pyarrow_dtype') or str(obs_df[col].dtype) == 'string[pyarrow]':
            obs_df[col] = obs_df[col].astype(str).astype(object)
        elif not pd.api.types.is_numeric_dtype(obs_df[col]) and not isinstance(obs_df[col].dtype, pd.CategoricalDtype):
            obs_df[col] = obs_df[col].astype(str).astype(object)
    new_index = obs_df.index.astype(str).astype(object)
    obs_df.index = new_index
    adata.obs = obs_df
    var_df = adata.var.copy()
    for col in var_df.columns:
        if hasattr(var_df[col].dtype, 'pyarrow_dtype') or str(var_df[col].dtype) == 'string[pyarrow]':
            var_df[col] = var_df[col].astype(str).astype(object)
        elif not pd.api.types.is_numeric_dtype(var_df[col]) and not isinstance(var_df[col].dtype, pd.CategoricalDtype):
            var_df[col] = var_df[col].astype(str).astype(object)
    new_var_index = var_df.index.astype(str).astype(object)
    var_df.index = new_var_index
    adata.var = var_df
    return adata

adata = fix_anndata_definitively(adata)
adata.write("adata_immuneFiltered.h5ad")
