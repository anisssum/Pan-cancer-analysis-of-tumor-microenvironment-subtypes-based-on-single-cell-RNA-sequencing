import os
import numpy as np
import pandas as pd
import scanpy as sc

# Paths
folder_path = 'Path to the folder atlas_dataset'
h5ad_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.h5ad')]
metadata = pd.read_csv('data_bulk.csv')
metadata_cells = metadata['cell_id'].tolist()

# Define genes of interest
genes_of_interest = ['PTPRC', 'CD4']

gene_stats = {
    'gene1_only': 0,
    'gene2_only': 0,
    'both_present': 0,
    'none_present': 0
}

def transform_values(matrix):
    non_zero_indices = matrix != 0
    matrix[non_zero_indices] = np.round(np.exp(matrix[non_zero_indices]))
    return matrix

# Process each .h5ad file
for h5ad_file in h5ad_files:
    adata = sc.read_h5ad(h5ad_file)
    dataset_name = os.path.basename(h5ad_file)
    print(f"Processing dataset: {dataset_name}")

    sample_ids = metadata['Sample'].unique()
    for sample_id in sample_ids:
        sample_indices = adata.obs.index[(adata.obs['Sample'] == sample_id) & (adata.obs.index.isin(metadata_cells))]
        if len(sample_indices) > 1:
            expr_matrix = adata[sample_indices].X.T 
            genes = adata.var_names
            cells = sample_indices

            expr_matrix = transform_values(expr_matrix.toarray())
            expr_df = pd.DataFrame(expr_matrix, index=genes, columns=cells)

            gene1_expr = expr_df.loc[genes_of_interest[0]] if genes_of_interest[0] in expr_df.index else pd.Series(0, index=cells)
            gene2_expr = expr_df.loc[genes_of_interest[1]] if genes_of_interest[1] in expr_df.index else pd.Series(0, index=cells)

            for cell in cells:
                gene1_value = gene1_expr[cell]
                gene2_value = gene2_expr[cell]

                if gene1_value > 0 and gene2_value > 0:
                    gene_stats['both_present'] += 1
                elif gene1_value > 0:
                    gene_stats['gene1_only'] += 1
                elif gene2_value > 0:
                    gene_stats['gene2_only'] += 1
                else:
                    gene_stats['none_present'] += 1


print('Gene Expression Statistics:')
print(gene_stats)
