######## Load all libraries 

#single cell libraries
import scanpy as sc
import anndata as ad


from tqdm import tqdm
import os
import pandas as pd
import warnings
import time
import argparse
warnings.simplefilter(action='ignore', category=FutureWarning)

## Load scripts with functions
import sys
sys.path.append("/gstore/project/elenakharitonova_melocars/downsampling_perturbseq/Pipeline_Development/utils")
import scProc
import rapids_modified_10000
import normalization_replogle



##########################################################################################################
########################### Part 1: Load Arguments and Data ##############################################
##########################################################################################################

## Load data
adata = sc.read_h5ad("/gstore/project/elenakharitonova_melocars/downsampling_perturbseq/Pipeline_Development/data/Sublib1_Sublib2_Sublib4_DS000017114_DS000016652_remove_pos_cont_counts_obs_var.h5ad")
print("loaded concat data")

## Create a numeric column in var, so that cudf series can be created
adata.var['num_id'] = range(len(adata.var))
adata.var['Symbol'] = adata.var.index

###################################################################################################
########################### Part 2: Preprocess Data and Obtain 2000 most variable genes ###########################################
###################################################################################################

## Create object to undergo preprocessing 
RP = rapids_modified_10000.RapidsSingleCellPipeline(adata)

## Preprocess data
RP.proc(10000, norm=True, scale=True, regress=False, embedding=True,
                 n_components=50, n_neighbors=10, knn_n_pcs=30, batch_key=None, filtered= False)

## Add back info about gene ids
RP.adata.var = RP.adata_.var[RP.adata_.var["num_id"].isin(RP.adata.var.index.to_list())]

## Obtain list of top 2000 most variable genes
HVG=RP.adata.var.index.tolist()

## Save list
HVG_file_out = "/gstore/scratch/u/kharitoe/sublib1_2_4/data/nopos_controls_sublib1_2_4_HVG.csv"
RP.adata.var.to_csv(HVG_file_out, index = False)

## Save Sample Information 
sample_file_out = "/gstore/scratch/u/kharitoe/sublib1_2_4/data/nopos_controls_sublib1_2_4sample_info.csv"
adata.obs.to_csv(sample_file_out)

print("Done Calculating HVG")
