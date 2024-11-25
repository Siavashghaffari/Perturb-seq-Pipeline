import scanpy as sc
import numpy as np
import scipy
import pandas as pd
from scipy.io import mmwrite
import argparse
import os


file_load =  "/gstore/scratch/u/ghaffars/glmGamPoi/sublib1_bdev/data/remove_pos_cont_counts_obs_var.h5ad"
data = sc.read_h5ad(file_load)

## Create NTC vs not label
data.obs["label"] = ['aNTC' if gene_symbol == "NTC" else 'zAll' for gene_symbol in data.obs['gene_symbol']]

data.obs['gem'] = data.obs['NGS_ID'].astype(str) + '-' + data.obs['10Xrun'].astype(str)


## Meta Data for pyGamPoi Model
meta = pd.DataFrame({'label':data.obs["label"],
                   'cell':data.obs.index,
                   'total_counts':data.obs["total_counts"],
                   'fct_counts_mt':data.obs["pct_counts_mt"]/100,
                   'batchid':data.obs["gem"],
                   'guidename':data.obs["DemuxAssignment_crispr"],
                   'targetname':data.obs["gene_symbol"]})


## Change guidename to not have commas to not confuse bash later
meta["guidename"] = meta["guidename"].str.replace(",", ":")


## List of targets
targets = data.obs["gene_symbol"].unique()

## Create Folder To Save by Target
by_target_folder =  "/gstore/scratch/u/ghaffars/glmGamPoi/sublib1_bdev/data/by_target/"

## Ensure folder exists
if not os.path.exists(by_target_folder):
    os.makedirs(by_target_folder)

## Save info about column names of top 3000 genes
print("length target is " + str(len(targets)))

## Iterate over target and save    
for i in range(len(targets)): 
    if i % 100 == 0:
            print(i)
    target = targets[i]
    data_target = data[data.obs["gene_symbol"].isin([target,"NTC"])].copy()
    data_target.X = data_target.X.astype("float32")
    data_target.layers['counts'] = data_target.layers['counts'].astype("float32")
    ## Save as .h5ad
    data_target.write(by_target_folder+target+".h5ad")
    

print("filtered and saved by target")
print("all done")