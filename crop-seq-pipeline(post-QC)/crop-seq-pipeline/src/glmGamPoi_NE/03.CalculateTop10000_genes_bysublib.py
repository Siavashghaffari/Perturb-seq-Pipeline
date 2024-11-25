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

## Load arguments
parser = argparse.ArgumentParser(description='Script to downsample depth and samples')
parser.add_argument('-s','--sublib', required =True,
                    help='Which sublib to analyze (h5ad file)')
args = parser.parse_args()
print(args)

## Specify input file depending on sublib
if args.sublib == "sublib1":
    file = "/gstore/project/elenakharitonova_melocars/downsampling_perturbseq/Pipeline_Development/data/Sublib1_DS000017114_remove_pos_cont_counts_obs_var.h5ad"
if args.sublib == "sublib2":
    file = "/gstore/scratch/u/kharitoe/Sublib2_Sublib4_DS000016652_no_pos_cont_downsample/data/downsample_data/downsampledepth1.0_downsamplesample1.0_seed1_rep0_sublib2.h5ad"
if args.sublib == "sublib4":
    file = "/gstore/scratch/u/kharitoe/Sublib2_Sublib4_DS000016652_no_pos_cont_downsample/data/downsample_data/downsampledepth1.0_downsamplesample1.0_seed1_rep0_sublib4.h5ad"

## Load data
adata = sc.read_h5ad(file)
print("loaded " + arg.sublib + " data")

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
HVG_file_out = "/gstore/scratch/u/kharitoe/sublib1_2_4/data/nopos_controls_" + args.sublib + "_top10000_HVG.csv"
RP.adata.var.to_csv(HVG_file_out, index = False)

## Save Sample Information 
sample_file_out = "/gstore/scratch/u/kharitoe/sublib1_2_4/data/nopos_controls_" + args.subib + "_sample_info.csv"
adata.obs.to_csv(sample_file_out)

print("Done Preprocessing")