#!/usr/bin/env python
# srun --time=72:00:00 --mem=500G --cpus-per-task=8 --qos=achesnokova --pty bash
import scanpy as sc
import pandas as pd
import numpy as np
import gc
import os
import harmonypy as hm
import celltypist
from celltypist import models
import anndata

# =========================================================
# SETTINGS
# =========================================================

pd.options.future.infer_string = False
anndata.settings.allow_write_nullable_strings = True

sc.settings.n_jobs = 8

OUTPUT_PREFIX = "/mnt/tank/scratch/achesnokova/flamingo/korean_tisch_3ca_data/analysis_results_harmony2/"
os.makedirs(OUTPUT_PREFIX, exist_ok=True)

# =========================================================
# LOAD DATA
# =========================================================

adata = sc.read_h5ad(
    "/mnt/tank/scratch/achesnokova/flamingo/korean_tisch_3ca_data/merged_all_464.h5ad"
)
# =========================================================
# MALIGNANT PSEUDOBULK
# =========================================================
from scipy.sparse import issparse

malignant = adata[
    adata.obs["celltype"] == "Malignant"
]

sample_ids = malignant.obs["sample_id"].astype(str)
unique_samples = sample_ids.unique()

pb_dict = {}

for sample in unique_samples:

    idx = np.where(sample_ids == sample)[0]

    X = malignant.X[idx]

    if issparse(X):
        pb_dict[sample] = np.asarray(
            X.sum(axis=0)
        ).ravel()
    else:
        pb_dict[sample] = np.array(
            X.sum(axis=0)
        ).ravel()

pb = pd.DataFrame(
    pb_dict,
    index=malignant.var_names
)

pb.to_csv(
    f"{OUTPUT_PREFIX}pseudobulk_malignant.csv"
)
# =========================================================
# INITIAL COMPOSITION MATRIX
# =========================================================

all_samples = adata.obs["sample_id"].astype(str).unique()

malignant_counts = (
    malignant.obs["sample_id"]
    .astype(str)
    .value_counts()
    .reindex(all_samples)
    .fillna(0)
    .astype(int)
)

composition_matrix = pd.DataFrame(
    [malignant_counts.values],
    index=["Malignant"],
    columns=all_samples
)

del malignant
del pb
gc.collect()

# =========================================================
# REMOVE MALIGNANT CELLS
# =========================================================
adata = adata[
    adata.obs["celltype"] != "Malignant"
].copy()

print(f"Cells: {adata.shape[0]}, Genes: {adata.shape[1]}")

gc.collect()

# =========================================================
# GENE FILTERING
# =========================================================

genes = pd.read_csv(
    "/mnt/tank/scratch/achesnokova/flamingo/korean_tisch_3ca_data/genes_464.txt",
    sep="\t"
)["x"].tolist()

genes = [
    g for g in genes
    if not g.startswith("MT-")
]

genes = [
    g for g in genes
    if g in adata.var_names
]

adata = adata[:, genes].copy()

sc.pp.filter_genes(
    adata,
    min_cells=10
)

print(f"After filtering: {adata.shape}")

# =========================================================
# NORMALIZATION
# =========================================================

sc.pp.normalize_total(
    adata,
    target_sum=1e4
)

sc.pp.log1p(adata)

# =========================================================
# CREATE STABLE BATCHES
# =========================================================

sample_sizes = adata.obs["sample_id"].value_counts()

small_samples = sample_sizes[
    sample_sizes < 50
].index

adata.obs["integration_batch"] = (
    adata.obs["sample_id"]
    .astype(str)
)

adata.obs.loc[adata.obs["sample_id"].isin(small_samples), "integration_batch"] = "SMALL_BATCH"

adata.obs["integration_batch"] = (
    adata.obs["integration_batch"]
    .astype("category")
)

print(f"Small samples merged: {len(small_samples)}")
print(small_samples)
# =========================================================
# HVGs
# =========================================================

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=5000,
    flavor="seurat_v3",
    batch_key="integration_batch"
)

adata = adata[:,adata.var.highly_variable].copy()

print(f"HVG matrix: {adata.shape}")

# =========================================================
# PCA
# =========================================================

sc.tl.pca(
    adata,
    n_comps=100,
    svd_solver="randomized",
    zero_center=False
)

# =========================================================
# HARMONY2
# =========================================================

ho = hm.run_harmony(
    adata.obsm["X_pca"],
    adata.obs,
    "integration_batch",
    theta=4.0,
    sigma=0.1,
    max_iter_harmony=10,
    verbose=True
)

adata.obsm["X_pca_harmony"] = ho.Z_corr
del ho
gc.collect()
# =========================================================
# NEIGHBORS + UMAP
# =========================================================

sc.pp.neighbors(
    adata,
    use_rep="X_pca_harmony",
    n_neighbors=30)

sc.tl.umap(
    adata,
    min_dist=0.1
)

sc.tl.leiden(
    adata,
    resolution=1,
    flavor="igraph",
    n_iterations=2
)
# =========================================================
# CELLTYPIST
# =========================================================

models.download_models(
    model="Immune_All_Low.pkl"
)

pred = celltypist.annotate(
    adata,
    model="Immune_All_Low.pkl",
    majority_voting=True,
    over_clustering="leiden"
)

adata.obs["celltypist_predicted_labels"] = pred.predicted_labels["majority_voting"]

adata.obs["celltypist_conf_score"] = pred.probability_matrix.max(axis=1).values

# =========================================================
# CELL COMPOSITION MATRIX
# =========================================================

ct_table = pd.crosstab(
    adata.obs["celltypist_predicted_labels"],
    adata.obs["sample_id"]
)

ct_table = ct_table.reindex(
    columns=all_samples,
    fill_value=0
)

composition_matrix = pd.concat(
    [composition_matrix, ct_table],
    axis=0
)

composition_matrix = composition_matrix.astype(int)

composition_matrix.to_csv(
    f"{OUTPUT_PREFIX}cell_composition_matrix.csv"
)
# =========================================================
# FIX STRINGS FOR WRITING
# =========================================================

def fix_anndata_for_writing(adata):
    obs_dict = {}
    for col in adata.obs.columns:
        series = adata.obs[col]
        if isinstance(
            series.dtype,
            pd.CategoricalDtype
        ):
            cats = (
                series.cat.categories
                .astype(str)
                .astype(object)
            )
            obs_dict[col] = (
                pd.Categorical.from_codes(
                    series.cat.codes,
                    categories=cats,
                    ordered=series.cat.ordered
                )
            )
        elif pd.api.types.is_numeric_dtype(series):
            obs_dict[col] = series.values
        else:
            obs_dict[col] = (
                series.astype(str)
                .astype(object)
                .values
            )
    adata.obs = pd.DataFrame(
        obs_dict,
        index=adata.obs.index.astype(str)
    )
    var_dict = {}
    for col in adata.var.columns:
        series = adata.var[col]
        if isinstance(
            series.dtype,
            pd.CategoricalDtype
        ):
            cats = (
                series.cat.categories
                .astype(str)
                .astype(object)
            )
            var_dict[col] = (
                pd.Categorical.from_codes(
                    series.cat.codes,
                    categories=cats,
                    ordered=series.cat.ordered
                )
            )
        elif pd.api.types.is_numeric_dtype(series):

            var_dict[col] = series.values
        else:
            var_dict[col] = (
                series.astype(str)
                .astype(object)
                .values
            )
    adata.var = pd.DataFrame(
        var_dict,
        index=adata.var.index.astype(str)
    )
    return adata

adata = fix_anndata_for_writing(adata)

if "highly_variable" in adata.var.columns:
    adata.var["highly_variable"] = (
        adata.var["highly_variable"]
        .astype(bool)
    )

# =========================================================
# SAVE
# =========================================================

final_df = adata.obs[
    [
        "celltypist_predicted_labels",
        "celltypist_conf_score"
    ]
].copy()

final_df[
    ["UMAP1", "UMAP2"]
] = adata.obsm["X_umap"]

final_df.to_csv(
    f"{OUTPUT_PREFIX}final_annotations.tsv",
    sep="\t"
)

adata.write(
    f"{OUTPUT_PREFIX}annotated_adata.h5ad"
)
# =========================================================
# PLOTS
# =========================================================
import matplotlib.pyplot as plt

sc.settings.set_figure_params(
    dpi=150,
    dpi_save=300,
    fontsize=10
)

# ---------- Cell types ----------
plt.figure(figsize=(24, 12))
sc.pl.umap(
    adata,
    color="celltypist_predicted_labels",
    size=1,
    legend_loc='right margin',
    legend_fontsize=8,
    frameon=False,
    title='CellTypist labels',
    show=False
)
plt.savefig(
    f"{OUTPUT_PREFIX}umap_harmony_celltype.png",
    bbox_inches='tight'
)
plt.close()

plt.figure(figsize=(24, 12))
sc.pl.umap(
    adata,
    color="Cancer type",
    size=1,
    legend_loc='right margin',
    legend_fontsize=8,
    frameon=False,
    title='Cancer type',
    show=False
)
plt.savefig(
    f"{OUTPUT_PREFIX}umap_harmony_cancer_type.png",
    bbox_inches='tight'
)
plt.close()

# =========================================================
# 3. CROSS ANNOTATION MATRIX
# =========================================================

OLD_COL = "celltype"
NEW_COL = "celltypist_predicted_labels"

ct_old = adata.obs[OLD_COL].astype(str)
ct_new = adata.obs[NEW_COL].astype(str)

# cross table
cross = pd.crosstab(ct_old, ct_new)
out3_csv = os.path.join(OUTPUT_PREFIX, "annotation_comparison_counts.csv")
cross.to_csv(out3_csv)
