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
      




