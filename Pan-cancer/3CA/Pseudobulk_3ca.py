import os
import pandas as pd
from glob import glob

root_dir = "/home/anya/flamingo/3CA/"

all_dfs = []

for folder in os.listdir(root_dir):

    folder_path = os.path.join(root_dir, folder)
    if not os.path.isdir(folder_path):
        continue
    print(folder)
    matrix_files = glob(os.path.join(folder_path, "sum_matrix_results*.csv"))

    for file_path in matrix_files:
        df = pd.read_csv(file_path, index_col=0)
        if not df.index.is_unique:
            print("dublicates in ", folder)
            df = df.groupby(df.index).sum()
        
        all_dfs.append(df)

merged_matrix = pd.concat(all_dfs, axis=1, join="outer")

merged_matrix.to_csv("Merged_Gene_Expression_Matrix_2.csv")
