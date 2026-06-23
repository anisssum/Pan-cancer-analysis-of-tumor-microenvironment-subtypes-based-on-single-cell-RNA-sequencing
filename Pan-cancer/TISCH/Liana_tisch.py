import os
import scanpy as sc
import pandas as pd
import liana as li
from pathlib import Path

data_dir = Path("tish_data_h5")
samples_keep = pd.read_csv("/home/anya/flamingo/samples_464.txt", header=None)[0].tolist()

meta = pd.read_csv("Merged_Cells_Final_Annotated_TISCH.txt", sep="\t")

output_dir = "liana_results"
os.makedirs(output_dir, exist_ok=True)

meta = meta[meta['Sample_fin'].isin(samples_keep)]

cells_types_list = [
    "B cell", "CD4", "CD8", "Dendritic cell", "Endothelial",
    "Fibroblast", "Macrophage", "Malignant", "Mast",
    "NK cell", "Plasma cell", "Treg", "Unknown_T_cells"
]

meta = meta[meta["Celltype_fin"].isin(cells_types_list)]

celltype_map = dict(zip(meta["Cell"], meta["Celltype_fin"]))
sample_map   = dict(zip(meta["Cell"], meta["Sample_fin"]))

for h5_file in data_dir.glob("*.h5"):
    try:
        adata = sc.read_10x_h5(h5_file, gex_only=False)

        dataset_results = []

        sample_prefix = Path(h5_file).stem.split('_expression')[0]

        current_sample_cells = meta[meta["Dataset"].str.startswith(sample_prefix)]

        adata = adata[current_sample_cells["Cell"], :].copy()

        adata.obs["Celltype_fin"] = adata.obs.index.map(celltype_map)
        adata.obs["Sample"] = adata.obs.index.map(sample_map)

        adata.layers["log"] = adata.X.copy()

        for sample in adata.obs["Sample"].unique():
            try:
                adata_sample = adata[adata.obs["Sample"] == sample].copy()
                adata_sample.layers["log"] = adata_sample.X.copy()
                li.mt.rank_aggregate(
                    adata=adata_sample,
                    groupby="Celltype_fin",
                    expr_prop=0.05,
                    min_cells=3,
                    resource_name="consensus",
                    use_raw=False,
                    layer="log",
                    verbose=True
                )

                res = adata_sample.uns["liana_res"]
                res["sample"] = sample
                res["dataset"] = sample_prefix
                dataset_results.append(res)

            except Exception as e:
                print(f" {sample} in {sample_prefix}: {e}")
                continue
              
        if dataset_results:
            liana_all = pd.concat(dataset_results)
            output_path = os.path.join(output_dir, f"liana_results_{sample_prefix}.csv")
            liana_all.to_csv(output_path, index=False)

    except Exception as e:
        print(f"{h5_file.name}: {e}")
        continue
