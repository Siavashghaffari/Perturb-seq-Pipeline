import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import os
import pickle

with open("/gstore/scratch/u/ghaffars/Dataset/Norm_control/Sublib4/HVG.pkl", "rb") as fp:   # Unpickling
    HVG = pickle.load(fp)

folder_path= "/gstore/project/crc_recursion_gw/DLD1_Sublib4/DS000016763/test_norm/"
file_names=[]
file_paths=[]
for root,dirs, files in os.walk(folder_path):
    for file in files:
        file_path = os.path.join(root,file)
        file_names.append(file)
        file_paths.append(file_path)
        
adatas=[]
for file in file_paths:
    AD = sc.read(file)
    print(file)
    AD.X.data = np.nan_to_num(AD.X.data, copy=False)
    AD.X.data=AD.X.data.astype("float32")
    AD = AD[:,AD.var.index.isin(HVG)].copy()
                
    adatas.append(AD)
    

adata = ad.concat(adatas,join='outer')
adata.var["Symbol"] = adata.var.index.copy()

adata.write("/gstore/scratch/u/ghaffars/Dataset/Norm_control/Sublib4/concat_1.h5ad")