import os
import scanpy as sc
import pandas as pd
import liana as li
import numpy as np
import matplotlib.pyplot as plt

# Paths and configurations
input_dir = "/home/anya/flamingo/Korean/atlas_dataset/"
output_dir = "/home/anya/flamingo/Korean/liana_results/"

sample_list = pd.read_csv("/home/anya/flamingo/samples_464.txt", header=None)[0].tolist()
merged_annot = pd.read_csv("/home/anya/flamingo/Korean/Merged_Cells_Final_Annotated_Korean.txt", sep="\t")
celltype_map = dict(zip(merged_annot["cell_id"], merged_annot['CellTypeFinAnnot']))

# Process each h5ad file
for filename in os.listdir(input_dir):
    if filename.endswith(".h5ad"):
        dataset_results = []
        dataset_name = filename.replace(".h5ad", "")
        print(f"\n--- Processing {filename} ---")

        try:
            # Load the dataset
            adata = sc.read_h5ad(os.path.join(input_dir, filename))

            # Filter samples
            adata = adata[adata.obs["Sample"].isin(sample_list)].copy()
            if adata.n_obs == 0:
                print(f"  No matching samples in {dataset_name}. Skipping.")
                continue

            print(f"  Retained {adata.n_obs} cells after filtering.")

            adata.obs["CellTypeFinAnnot"] = adata.obs.index.map(celltype_map)

            # Run LIANA for each sample separately
            for sample in adata.obs["Sample"].unique():
                print(f"  Processing sample {sample}...")

                try:
                    adata_sample = adata[adata.obs["Sample"] == sample].copy()

                    if adata_sample.obs["CellTypeFinAnnot"].nunique() < 2:
                        print(f"    Skipping: only one cell type.")
                        continue

                    # ensure layer exists
                    adata_sample.layers["log"] = adata_sample.X.copy()

                    # LIANA run
                    li.mt.rank_aggregate(
                        adata=adata_sample,
                        groupby="CellTypeFinAnnot",
                        expr_prop=0.05,
                        min_cells=3,
                        resource_name="consensus",
                        use_raw=False,
                        layer="log",
                        verbose=True
                    )

                    res = adata_sample.uns["liana_res"]
                    res["sample"] = sample
                    res["dataset"] = dataset_name
                    dataset_results.append(res)

                except Exception as sample_error:
                    print(f"    ERROR in sample {sample}: {sample_error}")
                    print("    → Skipping this sample and continuing to the next one.")
                    continue  # ← skip only this sample, continue dataset

            # Save combined results for dataset
            if dataset_results:
                liana_all = pd.concat(dataset_results)
                output_path = os.path.join(output_dir, f"liana_results_{dataset_name}.csv")
                liana_all.to_csv(output_path, index=False)
                print(f"Saved results for {dataset_name} to {output_path}")

        except Exception as dataset_error:
            print(f"DATASET ERROR for {filename}: {dataset_error}")
            print("→ Skipping entire dataset.")
            continue
