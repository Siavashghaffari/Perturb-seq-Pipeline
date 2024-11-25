## Import Libraries

from __future__ import division 

#single cell libraries
import scanpy as sc
import anndata
#import scvi

#general
import pandas as pd
import sys, os
import numpy as np
import itertools
from glob import glob
import warnings
from tqdm import tqdm
warnings.simplefilter(action='ignore', category=FutureWarning)
from collections import OrderedDict 
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

#DatasetDB
import pydsdb
from pydsdb import get_datasets

# Helper Scripts
import Templates.tools.scProc as proc

def load_data(adata,  alt_experiments, Index=False):
    """
    This function just cleans the data
    """
    # Change index in var to symbol
    try: 
        adata.var.index = adata.var['Symbol'].apply(lambda x:x.split('_')[1])
    except IndexError:
        if Index:
            adata.var.index = adata.var['Symbol']
        
    adata.var_names_make_unique()
    if len(alt_experiments)==1:
        adata.obs.columns = ["Sample", "Barcode", "DemuxType_hashing", "DemuxAssignment_hashing"]
    return adata


def Summerize_counts (adata, groupby):
    """
    The function summerizes the counts
    """
    cnt = proc.count_features(adata.obs, [groupby]).set_index([groupby][0])
    return cnt


def upper_lower_limits (adata, groupby):
    # Define varibale which come handy later
    variables = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo']
    # Here the limits that define bottom outliers
    LOWER = pd.concat([pd.DataFrame.from_dict(
        proc.calc_lower_limit_per_group(adata, groupby, el), orient='index').rename(
        columns={0:el}) for el in variables], axis=1)
    #Here the limits that define upper outliers
    UPPER = pd.concat([pd.DataFrame.from_dict(
        proc.calc_upper_limit_per_group(adata, groupby, el), orient='index').rename(
            columns={0:el}) for el in variables], axis=1)
    return LOWER,UPPER

    
def melt(df, groupby, variables, step):
    df = pd.melt(df, id_vars=groupby, value_vars=variables)
    df['step'] = step
    return df


def aggregation (adata, groupby, agg, step):
    # Define varibale which come handy later
    variables = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo']
    if step=="before":
        aggr = melt(adata.obs, groupby, variables, step)
    elif step=="after":
        aggr = pd.concat([agg, melt(adata.obs, groupby, variables, step)])
    return aggr


def doublet_counter (adata):
    total = adata.obs.shape[0]
    singlet_hashing = adata[adata.obs['DemuxType_hashing']=='singlet'].shape[0]
    if 'DemuxType_crispr' in adata.obs.columns:
        singlet_crispr = adata[adata.obs['DemuxType_crispr']=='singlet'].shape[0]
        singlet_both = adata[(adata.obs['DemuxType_crispr']=='singlet') & (adata.obs['DemuxType_hashing']=='singlet')].shape[0]
        d = {'Cell counts': [total, singlet_hashing, singlet_crispr, singlet_both]}
        df = pd.DataFrame(data=d, index=["total Cells", "cells after keeping singlets by hashing", 
    "cells after keeping singlets by gRNA", "cells after keeping singlets by both hashing and gRNA"])
    else:
        d = {'Cell counts': [total, singlet_hashing]}
        df = pd.DataFrame(data=d, index=["total Cells", "cells after keeping singlets by hashing"])
    return df


def pipeline_has_gpu():
    try:
        import cupy as cp
        return True
    except ImportError:
        return False
    
    
def Mapping_reads_info (Mapping_rate):
    
    # Create a dictionary to store mapping rates of different libraries
    extraced_data={}
    
    for key,value in Mapping_rate.items():
        #split the value by the newline character
        parts = value.split("\n")
        # remove empty str from the list
        parts =[part for part in parts if part !='']
        # create an empty dictionary to store extracted values
        values_dict={}
        for part in parts:
            try:
                label,count = part.split(":")
            except ValueError:
                Label = part
            values_dict[label]=count
        values_dict["Experiment"]=Label.split(".")[1]
        # Create adtafarme from values dictionary
        df = pd.DataFrame(values_dict, index=[key])
        #Store the daatafrem in the datafarme dictionary
        extraced_data[key] = df

    #display tables
    combined_df = pd.concat(extraced_data.values())
    agg = combined_df[["Total number of reads","Number of reads with valid cell and feature barcodes","Number of valid UMIs (with matching cell and feature barcodes)"]]
    agg.columns=["Total Reads", "Mapped Reads", "UMIs"]
    agg["Total Reads"]=agg["Total Reads"].apply(lambda x:int(x))
    agg["Mapped Reads"]=agg["Mapped Reads"].apply(lambda x:int(x.split("(")[0]))
    agg["UMIs"]=agg["UMIs"].apply(lambda x:int(x))
    #Make a dataframe for percentage data as well
    agg_percentage = agg.copy()
    agg_percentage["Total Reads"]=agg["Total Reads"]/agg["Total Reads"]*100
    agg_percentage["Mapped Reads"]=agg["Mapped Reads"]/agg["Total Reads"]*100
    agg_percentage["UMIs"]=agg["UMIs"]/agg["Total Reads"]*100
    return combined_df,agg,agg_percentage

def Ven_maker(adata, Title):
    crispr = set(adata.obs[adata.obs['DemuxType_crispr']=='singlet'].index.to_list())
    hashing = set(adata.obs[adata.obs['DemuxType_hashing']=='singlet'].index.to_list())
    total = len(hashing.union(crispr))
    venn2([hashing, crispr],set_labels = ('singlets by hashing', 'singlets by crispr'),
     subset_label_formatter=lambda x: str(x) + "\n(" + f"{(x/total):1.0%}" + ")")
    plt.title(Title)
    plt.show()