#single cell libraries
import scanpy as sc
import anndata
#import scvelo as scv
#import scvi
from harmony import harmonize
#import muon as mu
# Import a module with ATAC-seq-related functions
#from muon import atac as ac

#general
import pandas as pd
import sys, os
import numpy as np
import scipy
from scipy import stats
import itertools
from glob import glob
import warnings
from tqdm import tqdm
warnings.simplefilter(action='ignore', category=FutureWarning)
from collections import OrderedDict 
from scipy.stats import zscore
import joypy
from scipy.sparse import csr_matrix


#plotting
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib import cm

### Hashing functions

def normalize_clr(col):
    """ The centered log-ratio (clr) transformation uses the geometric mean of the sample vector as the reference """

    val_array = np.array(col.values) + 1
    col_gmean = stats.mstats.gmean(val_array)

    col_norm = [np.log(i / col_gmean) for i in val_array]
    
    return(pd.Series(data=col_norm, index=col.index))


def norm_matrix(HTO_mtx_raw):
    HTO_mtx_norm = HTO_mtx_raw.to_df().apply(func=normalize_clr, axis=1)
    HTO_mtx_norm_adata = anndata.AnnData(HTO_mtx_norm)
    return HTO_mtx_norm_adata


def procBC(adata, n_pcs=20, n_neighbors=30, **kwargs):
    print("applying center log ratio")
    adata.X = norm_matrix(adata).X
    print("computing pca")
    sc.tl.pca(adata, svd_solver='arpack')
    print("computing neighbors")
    sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=n_neighbors, **kwargs)
    print("computing umap")
    sc.tl.umap(adata, random_state=42)
    return adata
    

def meltdf(adata, features):
    """features is a list"""
    df = adata.to_df()[features]
    df = pd.melt(df, value_vars=features)
    return df


def plot_ridge(df, title, **kwargs):
    joypy.joyplot(
    data=df,
    by='BC',
    column='value',
    colormap=cm.tab20,
    title=title,
    **kwargs
    )
    
def count_features(adata, groupby):
    cnts = adata.obs.groupby(groupby)[groupby[0]].count().reset_index(name='cnt')
    # filter cells below 3
    cnts = cnts[cnts['cnt']>3]
    cnts['pct'] = cnts.groupby([groupby[0]])['cnt'].transform(lambda x : np.round(100*x/x.sum(), 1))
    return cnts


def read_demuxEM_map(f):
    df = pd.read_csv(f).groupby(['sample'])['index'].agg(list).reset_index()
    df.columns = ['sample', 'assignment']
    df['assignment'] = df['assignment'].apply(lambda x:sorted(x))
    df['assignment'] = df['assignment'].apply(lambda x:','.join(x))
    return {k:v for k,v in zip(df['assignment'], df['sample'])}

def fix_assignment(df):
    # This makes sure that the assignment is sorted
    df['assignment'] = df['assignment'].apply(lambda x:list(str(x).split(',')))
    df['assignment'] = df['assignment'].apply(lambda x:sorted(x))
    df['assignment'] = df['assignment'].apply(lambda x:','.join(x))
    return df['assignment']
      




# Set of functions to demultiplex barcode data based on margin

def max_value_per_row(arr):
    return np.max(arr, axis=1)

def dominant_col(arr):
    # np.argmax gives the indices of max values along the rows
    return np.argmax(arr, axis=1)

def second_largest_in_sparse_row(sparse_arr):
    if not isinstance(sparse_arr, csr_matrix):
        raise ValueError("Input sparse_arr must be a scipy.sparse.csr.csr_matrix")

    second_largest_values = [] 

    for i in range(sparse_arr.shape[0]): 
        row = sparse_arr.getrow(i).toarray()[0]
        row_without_zeros = row[row!=0]  # Remove zero entries
        if len(row_without_zeros) < 2:  # If less than two non-zero values
            second_largest_values.append(0)  # Or another value
            continue
        largest_indexes = np.argpartition(row_without_zeros, -2)[-2:]
        second_largest = row_without_zeros[largest_indexes[np.argsort(row_without_zeros[largest_indexes])[0]]]
        second_largest_values.append(second_largest)

    return np.array(second_largest_values)

def proc_array(adata, key, layer='counts'):
    # key is either crispr or hashing
    arr = adata.layers[layer]
    dominant_index = dominant_col(arr)
    dominant_column = np.array(adata.var_names)[dominant_index].flatten()
    dominant_value = max_value_per_row(arr).todense()
    second_value = second_largest_in_sparse_row(arr)
    # Build df
    df = pd.DataFrame(index=adata.obs.index)
    df[f'max_{key}'] = dominant_column
    df[f'max_{key}_value'] = dominant_value
    df[f'second_best_{key}_value'] = second_value
    return df

def assign_gRNA_identity(df, key, min_counts_cutoff, min_margin):
    df[f'L2FC_margin_{key}'] = np.log2(df[f'max_{key}_value']/df[f'second_best_{key}_value'])    
    df[f'marginType_{key}'] = 'doublet'
    df.loc[(df[f'max_{key}_value']<=min_counts_cutoff), f'marginType_{key}'] = 'unknown'
    df.loc[(df[f'L2FC_margin_{key}']>np.log2(min_margin)), f'marginType_{key}'] = 'singlet'
    return df

