import scanpy as sc
import numpy as np
import pandas as pd
import anndata
from pathlib import Path
from scipy.sparse import issparse

data_dir = Path("tish_data_h5")
samples_file = "samples_tisch.txt"
metadata_file = "merged_table.tsv"
output_file = "pseudobulk_malignant_3.csv"

with open(samples_file) as f:
    samples_keep = set(line.strip() for line in f if line.strip())

meta = pd.read_csv(metadata_file, sep="\t", low_memory=False)

malignant_labels = [
    "AC-like Malignant",
    "OC-like Malignant",
    "OPC-like Malignant",
    "NB-like Malignant",
    "MES-like Malignant",
    "NPC-like Malignant"
]

meta["Celltype_fin"] = np.where(
    meta["Celltype_(major-lineage)"].isin(malignant_labels),
    "Malignant",
    meta["Celltype_(major-lineage)"]
)

celltype_map = dict(zip(meta["Cell"], meta["Celltype_fin"]))
sample_map = dict(zip(meta["Cell"], meta["Sample"]))
patient_map = dict(zip(meta["Cell"], meta["Patient"]))
meta = 0

target_genes = ['PTPRC', 'CD3D', 'IL7R', 'CD8B', 'CD8A', 'NKG7', 'CD69', 'C1QA', 'FCN1', 'MZB1']

pseudobulk_dict = {}
total_cells = 0
total_cells_added = 0

for h5_file in data_dir.glob("*.h5"):
    dataset_prefix = h5_file.stem.split("_expression")[0]
    print(f"Обрабатываю {h5_file.name} с префиксом {dataset_prefix}...")

    adata = sc.read_10x_h5(h5_file, gex_only=False)
    X_raw = adata.X
    X_raw.data = np.expm1(X_raw.data)
    adata = anndata.AnnData(X=X_raw, obs=adata.obs.copy(), var=adata.var.copy())

    adata.obs["Celltype_fin"] = adata.obs.index.map(celltype_map)
    adata.obs["Sample"] = adata.obs.index.map(sample_map)
    adata.obs["Patient"] = adata.obs.index.map(patient_map) 

    adata.obs["Dataset_Sample"] = dataset_prefix + "_" + adata.obs["Patient"] + "_" + adata.obs["Sample"]

    adata = adata[adata.obs["Celltype_fin"] == "Malignant"]
    adata = adata[adata.obs["Dataset_Sample"].isin(samples_keep)]

    if adata.n_obs == 0:
        continue

    target_gene_indices = np.where(adata.var_names.isin(target_genes))[0]

    for sample_id in adata.obs["Dataset_Sample"].unique():
        mask = adata.obs["Dataset_Sample"] == sample_id
        sample_indices = np.where(mask)[0]

        if len(sample_indices) > 1:
            sample_data = adata[sample_indices].X

            if len(sample_indices) > 100:
                total_cells += len(sample_indices)

            if len(target_gene_indices) > 0:
                target_expr = sample_data[:, target_gene_indices]
                if issparse(target_expr):
                    zero_expression_cells = np.array(target_expr.sum(axis=1) == 0).flatten()
                else:
                    zero_expression_cells = np.all(target_expr == 0, axis=1)

                filtered_sample_data = sample_data[zero_expression_cells, :]

                if filtered_sample_data.shape[0] > 100:
                    if issparse(filtered_sample_data):
                        counts_sum = np.array(filtered_sample_data.sum(axis=0)).flatten()
                    else:
                        counts_sum = filtered_sample_data.sum(axis=0)

                    series_sample = pd.Series(counts_sum, index=adata.var_names)
                    pseudobulk_dict[sample_id] = series_sample
                    total_cells_added += filtered_sample_data.shape[0]

pseudobulk_df = pd.DataFrame.from_dict(pseudobulk_dict, orient="index")
pseudobulk_df.to_csv(output_file)

print(f"Total cells added: {total_cells_added}")
print(f"Total cells filtered out: {total_cells}")
