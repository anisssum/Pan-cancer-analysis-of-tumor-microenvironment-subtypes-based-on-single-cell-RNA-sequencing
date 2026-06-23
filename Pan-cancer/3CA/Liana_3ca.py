import os
import scanpy as sc
import pandas as pd
import liana as li
import numpy as np
from scipy import io
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix
import anndata as ad

input_dir = "/home/anya/flamingo/3CA/"
output_dir = "/home/anya/flamingo/3CA/liana_results/"

samples_id_final_to_keep = pd.read_csv("samples_464.txt", header=None)[0].tolist()

merged_annot = pd.read_csv(os.path.join(input_dir, "Merged_Cells_Final_Annotated_3CA_2.txt"), sep= "\t")

def load_single_group(path):
    mtx_files = [f for f in os.listdir(path) if f.endswith(".mtx")]
    mtx_file = mtx_files[0]
    mtx = io.mmread(os.path.join(path, mtx_file)).tocsr().T
    cells_files = [f for f in os.listdir(path) 
                   if f.startswith("Cells") and f.endswith(".csv")]
    if not cells_files:
        cells_files = [f for f in os.listdir(path) 
                      if "cell" in f.lower() and f.endswith(".csv")]
    if not cells_files:
        raise FileNotFoundError(f"No cells files found in {path}")
    cells_file = cells_files[0]
    cells_df = pd.read_csv(os.path.join(path, cells_file))
    if "cell_name" in cells_df.columns:
        cells = cells_df["cell_name"].tolist()
    elif "Cell" in cells_df.columns:
        cells = cells_df["Cell"].tolist()
    else:
        cells = cells_df.iloc[:, 0].tolist()
    genes_files = [f for f in os.listdir(path) 
                  if (f.startswith("Genes") or f.startswith("genes")) 
                  and (f.endswith(".csv") or f.endswith(".txt"))]
    if genes_files:
        genes_file = genes_files[0]
        if genes_file.endswith(".csv"):
            genes_df = pd.read_csv(os.path.join(path, genes_file), header=None)
        else:  # .txt
            genes_df = pd.read_csv(os.path.join(path, genes_file), header=None, sep='\t')
        genes = [str(g) for g in genes_df[0]]
    else:
        genes_df = pd.read_csv(os.path.join(os.path.dirname(path.rstrip('/')), "genes.txt"), header=None, names=["gene"])
        genes = [str(g) for g in genes_df["gene"]]
        
    adata = sc.AnnData(X=mtx)
    adata.var_names = genes
    adata.obs_names = cells
    print(f"Loaded group from {path}: {adata.shape[0]} cells, {adata.shape[1]} genes")
    return adata

# -----------------------------------------------------------
# Load dataset and merge with global annotation
# -----------------------------------------------------------
def load_dataset_from_folder(folder):
    print(f"Loading dataset from: {folder}")
    
    group_folders = [f for f in os.listdir(folder) 
                    if f.startswith("Group") and os.path.isdir(os.path.join(folder, f))]
    
    adata_list = []
    
    if group_folders:
        print(f"Found {len(group_folders)} group folders: {group_folders}")
        for g in sorted(group_folders):
            gpath = os.path.join(folder, g)
            try:
                group_adata = load_single_group(gpath)
                group_adata.obs["Group"] = g
                adata_list.append(group_adata)
            except Exception as e:
                print(f"Error loading group {g}: {str(e)}")
                continue
        
        if not adata_list:
            raise ValueError(f"No groups could be loaded from {folder}")
        
        adata = adata_list[0].concatenate(
            adata_list[1:], 
            batch_key="original_group", 
            batch_categories=[g.obs["Group"][0] for g in adata_list],
            index_unique=None  # Это ключевой параметр!
        )
    else:
        print(f"No group folders found. Trying to load directly from {folder}")
        try:
            adata = load_single_group(folder)
            adata.obs["Group"] = "Group1"
        except Exception as e:
            print(f"Error loading from {folder}: {str(e)}")
            raise
    
    dataset_name = os.path.basename(folder)
    print(f"Dataset name: {dataset_name}")
    
    dataset_samples = merged_annot[merged_annot["sample"].str.startswith(dataset_name + "_")]
    if len(dataset_samples) == 0:
        print(f"Warning: No annotations found for dataset {dataset_name}")
        adata.obs["CellTypeFinAnnot"] = ""
    else:
        annots = dataset_samples[["cell_name", "CellTypeFinAnnot", "sample"]]
        annots = annots.set_index("cell_name")
        adata.obs = adata.obs.join(annots, how="left")
        print(f"Annotations loaded for {len(dataset_samples)} cells")
    
    return adata

os.makedirs(output_dir, exist_ok=True)

datasets = [d for d in os.listdir(input_dir)
            if (d.startswith("Data") and 
                os.path.isdir(os.path.join(input_dir, d)) and 
                not os.path.exists(os.path.join(output_dir, f"liana_results_{d}.csv")))]


all_results = []
cell_types_list = ["B cell", "CD4", "CD8", "Dendritic cell", "Endothelial", "Fibroblast", "Macrophage", "Malignant", "Mast", "NK cell", "Plasma cell", "Treg", "Unknown_T_cells"]

def sum_duplicated_genes_sparse(adata):
    X = adata.X.tocsc()
    var_names = np.asarray(adata.var_names)
    unique_genes, inverse_idx = np.unique(var_names, return_inverse=True)
    rows = X.indices
    data = X.data
    col_ids = []
    for j in range(X.shape[1]):
        count = X.indptr[j+1] - X.indptr[j]
        if count > 0:
            col_ids.append(np.full(count, j, dtype=np.int32))
    col_ids = np.concatenate(col_ids)
    cols = inverse_idx[col_ids]
    X_new = coo_matrix(
        (data, (rows, cols)),
        shape=(X.shape[0], len(unique_genes))
    ).tocsr()
    new_adata = ad.AnnData(
        X=X_new,
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=unique_genes)
    )
    return new_adata

for dataset in datasets:
    dataset_path = os.path.join(input_dir, dataset)
    print(f"\nProcessing {dataset}...")
    try:
        adata = load_dataset_from_folder(dataset_path)
        adata = sum_duplicated_genes_sparse(adata)
        adata = adata[adata.obs['CellTypeFinAnnot'].isin(cell_types_list)].copy()
        if dataset == "Data_Hwang2022_Pancreas":
            sc.pp.log1p(adata)
        else:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        sc.pp.filter_genes(adata, min_cells=3)
        dataset_results = []
        unique_samples = adata.obs["sample"].unique()
        for sample in unique_samples:
            if sample not in samples_id_final_to_keep:
                continue
            adata_sample = adata[adata.obs["sample"] == sample].copy()
            if adata_sample.n_obs < 10:
                continue
            celltypes = adata_sample.obs["CellTypeFinAnnot"].unique()
            if len(celltypes) < 2:
                print(f"Skipping {sample}: only one cell type ({celltypes[0]})")
                continue
            print(f"Running LIANA for sample {sample} with {adata_sample.n_obs} cells...")
            li.mt.rank_aggregate(
                adata=adata_sample,
                groupby="CellTypeFinAnnot",
                expr_prop=0.05,
                min_cells=3,
                resource_name="consensus",
                use_raw=False,
                verbose=True
            )
            res = adata_sample.uns["liana_res"]
            res["sample"] = sample
            res["dataset"] = dataset
            dataset_results.append(res)
        if dataset_results:
            liana_all = pd.concat(dataset_results)
            unique_samples = liana_all["sample"].unique()
            missing_in_txt = [s for s in unique_samples if s not in samples_id_final_to_keep]
            extra_in_txt = [s for s in samples_id_final_to_keep if s not in unique_samples and s.startswith(dataset)]
            if missing_in_txt:
                print(f"Warning: Some samples not in samples_464.txt: {missing_in_txt}")
            if extra_in_txt:
                print(f"Info: Samples in samples_464.txt but not found in dataset: {extra_in_txt}")
            output_path = os.path.join(output_dir, f"liana_results_{dataset}.csv")
            liana_all.to_csv(output_path, index=False)
            print(f"Saved results for {dataset} to {output_path}")
            print(f"  Samples processed: {len(unique_samples)}")
        else:
            print(f"No samples from {dataset} were in samples_464.txt")
    except Exception as e:
        print(f"Error processing {dataset}: {str(e)}")
        import traceback
        traceback.print_exc()
        continue

if all_results:
    final_results = pd.concat(all_results)
    final_output_path = os.path.join(output_dir, "liana_results_all_datasets.csv")
    final_results.to_csv(final_output_path, index=False)
    final_samples = final_results["sample"].unique()
    print(f"\n=== FINAL SUMMARY ===")
    print(f"Saved combined results to {final_output_path}")
    print(f"Total interactions: {len(final_results)}")
    print(f"Unique sample in results: {len(final_samples)}")
    samples_in_both = [s for s in final_samples if s in samples_id_final_to_keep]
    samples_missing = [s for s in samples_id_final_to_keep if s not in final_samples]
    print(f"\nSamples successfully processed: {len(samples_in_both)}")
    if samples_missing:
        print(f"Samples from samples_464.txt that were NOT processed ({len(samples_missing)}):")
        for s in samples_missing[:20]:
            print(f"  - {s}")
        if len(samples_missing) > 20:
            print(f"  ... and {len(samples_missing) - 20} more")
    extra_samples = [s for s in final_samples if s not in samples_id_final_to_keep]
    if extra_samples:
        print(f"\nWarning: {len(extra_samples)} sample in results but not in samples_464.txt:")
        for s in extra_samples[:10]:
            print(f"  - {s}")
        if len(extra_samples) > 10:
            print(f"  ... and {len(extra_samples) - 10} more")
else:
    print("No results were generated!")

correspondence_df = pd.DataFrame({
    'sample_id_final_in_txt': samples_id_final_to_keep,
    'processed': [s in final_results['sample'].unique() if all_results else False for s in samples_id_final_to_keep]
})
correspondence_df.to_csv(os.path.join(output_dir, "samples_processing_status.csv"), index=False)
