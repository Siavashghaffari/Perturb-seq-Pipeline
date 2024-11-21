
#  This script Contains all helper functions necessary for the crop-seq pipeline
#    Modified by Siavash Ghaffari to be consistent with the crop-seq-pipeline


#single cell libraries
import scanpy as sc
import anndata as ad

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
import math
from scipy.sparse import csr_matrix


#plotting
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns

try:
    import cupy as cp
except ImportError:
    pass


def calc_qc_metrics(raw_adata, organism):
    mt = {'human':'MT-', 'mouse':'mt-'}
    ribo = {'human':("RPS","RPL"), 'mouse':("Rps","Rpl")}
    # save raw data
    raw_adata.layers["counts"] = raw_adata.X.copy()
    raw_adata.var['mt'] = raw_adata.var_names.str.startswith(mt[organism])
    raw_adata.var['ribo'] = raw_adata.var_names.str.startswith(ribo[organism])
    sc.pp.calculate_qc_metrics(raw_adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
    return raw_adata


def outliers(df):
    Q1=df.quantile(0.25)
    Q3=df.quantile(0.75)
    IQR=Q3-Q1
    upper = Q3+1.5*IQR
    lower = Q1-1.5*IQR
    return upper, lower


def calc_upper_limit_per_group(adata, groupby, x):
    d = {}

    for el in adata.obs[groupby].unique():
        d[el] = outliers(adata[adata.obs[groupby]==el].obs[x])[0]
    return d


def calc_lower_limit_per_group(adata, groupby, x):
    d = {}

    for el in adata.obs[groupby].unique():
        d[el] = outliers(adata[adata.obs[groupby]==el].obs[x])[1]
    return d


def filter_adata(adata, groupby, min_cells, min_genes, ribo=True, doublets=False, max_genes=None, pct_mt=None, total_counts=None):
    """
        If values are None, this triggers outlier detection for upper limit by group in groupby.
        min_cells and min_genes should be set by the user. 
        If ribo == True, lower and upper outliers are removed automatically
    
    """
    # Calculate outlier thresholds
    if max_genes == None:
        ge = calc_upper_limit_per_group(adata, groupby, 'n_genes_by_counts')
    if pct_mt == None:
        mt = calc_upper_limit_per_group(adata, groupby, 'pct_counts_mt')
    if total_counts == None:
        cn = calc_upper_limit_per_group(adata, groupby, 'total_counts')
    if ribo == True:
        ru = calc_upper_limit_per_group(adata, groupby, 'pct_counts_ribo')
        rl = calc_lower_limit_per_group(adata, groupby, 'pct_counts_ribo')
        
    print("starting with " + str(adata.X.shape[0]) + ' cells')
    
    sc.pp.filter_genes(adata, min_cells=int(min_cells))
    print("Now " + str(adata.X.shape[0]) + ' cells after filtering min cells ' + str(min_cells))
    
    sc.pp.filter_cells(adata, min_genes=int(min_genes))
    print("Now " + str(adata.X.shape[0]) + ' cells after filtering min genes ' + str(min_genes))
    
    if max_genes == None:
        print("Calculating upper outliers")
        for k in ge.keys():
            adata = adata[~((adata.obs[groupby]==k) & (adata.obs['n_genes_by_counts']>ge[k]))]
    else:
        adata = adata[(adata.obs['n_genes_by_counts']<=int(max_genes)).astype("bool")]
    print("Now " + str(adata.X.shape[0]) + ' cells after filtering max genes ')
    
    if pct_mt != None:
        adata = adata[adata.obs['pct_counts_mt'] < int(pct_mt)]
    else:
        print("Calculating upper outliers")
        for k in mt.keys():
            adata = adata[~((adata.obs[groupby]==k) & (adata.obs['pct_counts_mt']>mt[k]))]
    print("Now " + str(adata.X.shape[0]) + ' cells after filtering pct mt ')
    
    if total_counts!=None:
        adata = adata[adata.obs['total_counts']<=int(total_counts)]
    else:
        print("Calculating upper outliers")
        for k in cn.keys():
            adata = adata[~((adata.obs[groupby]==k) & (adata.obs['total_counts']>cn[k]))]
    print("Now " + str(adata.X.shape[0]) + ' cells after filtering total counts ')
    if ribo == True:
        for k in ru.keys():
            adata = adata[~((adata.obs[groupby]==k) & (adata.obs['pct_counts_ribo']>ru[k]))]
            adata = adata[~((adata.obs[groupby]==k) & (adata.obs['pct_counts_ribo']<rl[k]))]
        print("Now " + str(adata.X.shape[0]) + ' cells after filtering pct ribo ')
    if doublets == True:
        if 'DemuxType_hashing' in adata.obs.columns:
            adata = adata[adata.obs['DemuxType_hashing']=='singlet']
        elif 'DemuxType_crispr' in adata.obs.columns:
            adata = adata[adata.obs['DemuxType_crispr']=='singlet']
        print("Now " + str(adata.X.shape[0]) + ' cells after keeping singlets')
    return adata
    
    
def count_features(df, groupby):
    #df might adata.obs
    cnts = df.groupby(groupby)[groupby[0]].count().reset_index(name='cnt')
    # filter cells below 3
    cnts = cnts[cnts['cnt']>3]
    cnts['pct'] = cnts.groupby([groupby[0]])['cnt'].transform(lambda x : np.round(100*x/x.sum(), 1))
    return cnts

def pivot_features(cnts, index, columns, values):
    return pd.pivot_table(cnts, index=index, columns=columns, values=values)
    

def proc(adata, n_top_genes='auto', key=None, norm=True, scale=False, regress=False, embedding=True, n_pcs=30, n_neighbors=10, regress_cell_cycle=False, **hvg_kwargs):
    if norm == True:
        print('normalizing')
        adata.layers["counts"] = adata.X.copy() 
        sc.pp.normalize_total(adata, target_sum=10000)
        sc.pp.log1p(adata)
    if n_top_genes=='auto':
        print('selecting highly variable genes')
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, **hvg_kwargs)
        hvg=adata[:, adata.var.highly_variable].X.shape[1]
        print('Done selecting ' + str(hvg) + ' highly variable genes')
    else:
        print('selecting '+ str(n_top_genes) +' highly variable genes')
        sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=n_top_genes, min_mean=0.01, max_mean=5, min_disp=0.5, **hvg_kwargs)
    if regress==True:
        print("regressing")
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    if regress_cell_cycle==True:
        print("regressing cell cycle")
        sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
    if scale==True:
        print("scaling")
        sc.pp.scale(adata, max_value=10)
    if embedding==True:
        print("computing PCA")
        sc.tl.pca(adata)
        if key!=None:
            print('batch correcting')
            sc.external.pp.harmony_integrate(adata, basis='X_pca', adjusted_basis='X_pca_harmony', key=key, max_iter_harmony=50)
            sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_pcs=n_pcs, n_neighbors=30, random_state=42)
        else:
            sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=n_neighbors, random_state=42)
        print("computing UMAP")
        sc.tl.umap(adata, spread=1., min_dist=.5, random_state=11)
    return adata


def compute_cellCycleScores(adata, cellcyclefile):
    cell_cycle = pd.read_csv(cellcyclefile)
    s_genes = [x for x in cell_cycle[cell_cycle['phase']=='S']['symbol'].tolist() if x in adata.var_names]
    g2m_genes = [x for x in cell_cycle[cell_cycle['phase']=='G2/M']['symbol'].tolist() if x in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    return adata
    
    
def aggCnts(cnts_before, cnts_after):
    df = pd.concat([cnts_before[['cnt']], cnts_after[['cnt']]], axis=1)
    df.columns = ['before', 'after']
    #df = pd.melt(df.reset_index(), id_vars='ID', value_vars=['before', 'after'])
    df['pct left'] = np.round(100*df['after']/df['before'],1)
    return df


def plotViolin(data, variable, x, ax, **kwargs):
    sns.violinplot(data=data[data['variable']==variable], x=x, y='value', hue='step', inner="quartile", ax=ax, **kwargs)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
    for tick in ax.get_xticklabels():
        tick.set_rotation(90)
    sns.despine(ax=ax, top=True, right=True, left=False, bottom=False, offset=None, trim=False)

    
def plotQC(data, groupby, variable):
    figsize_x=len(data[groupby].unique())*1.1*2
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(figsize_x, figsize_x/3))
    plotViolin(data[data['step']=='before'], variable, groupby, ax1, palette='Blues')
    plotViolin(data[data['step']=='after'], variable, groupby, ax2, palette='Oranges')
    ax1.get_legend().remove()
    ax2.get_legend().remove()

    ax1.set_title('before', style='italic')
    f.suptitle(variable, y=1.1, x=0.5)
    ax2.set_title('after', style='italic')
    plt.tight_layout()
    
    
def plotQCSplit(data, groupby, variable):
    figsize_x=len(data[groupby].unique())*1.1
    f, ax1 = plt.subplots(1, figsize=(figsize_x, figsize_x/2))
    plotViolin(data, variable, groupby, ax1, split=True)
    ax1.set_title(variable, style='italic')
    plt.tight_layout()
    
    
def plot_cell_left(cnts_before, cnts_after, groupby):
    figsize_x=len(cnts_before)
    f, ax = plt.subplots(figsize=(figsize_x, 3))

    sns.barplot(data = aggCnts(cnts_before, cnts_after).reset_index(), y='pct left', x=groupby, color='tab:blue')
    for tick in ax.get_xticklabels():
        tick.set_rotation(90)
    ax.set_title('Percent cells left after QC filtering')
    
    
def getMean_genes(adata, category):
    res = pd.DataFrame(columns=adata.var_names, index=adata.obs[category].cat.categories)
    for clust in adata.obs[category].cat.categories: 
        res.loc[clust] = adata[adata.obs[category].isin([clust]),:].X.mean(0)
    return res


def dominant_col(df):
    return df.columns[np.argmax(df.values, axis=1)]


def plotbar(df, color, ax, **kwargs):
    df.plot(kind='bar', stacked=True, color=color, ax=ax, **kwargs)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

    

def Calculate_outlier_thresholds(adata, groupby, ribo=True, max_genes=None, pct_mt=None, total_counts=None): 
    # Calculate outlier thresholds
    if max_genes == None:
        ge = calc_upper_limit_per_group(adata, groupby, 'n_genes_by_counts')
    else:
        ge = None
    if pct_mt == None:
        mt = calc_upper_limit_per_group(adata, groupby, 'pct_counts_mt')
    else:
        mt = None
    if total_counts == None:
        cn = calc_upper_limit_per_group(adata, groupby, 'total_counts')
    else:
        cn = None
    if ribo == True:
        ru = calc_upper_limit_per_group(adata, groupby, 'pct_counts_ribo')
        rl = calc_lower_limit_per_group(adata, groupby, 'pct_counts_ribo')
    else:
        ru,rl = None, None
    return ge, mt, cn, ru, rl
    
    

def filter_adata_gpu(adata, groupby, min_cells, min_genes, ge, mt, cn, ru, rl, ribo=True, doublets=False, max_genes=None, pct_mt=None, total_counts=None):
    """
        If values are None, this triggers outlier detection for upper limit by group in groupby.
        min_cells and min_genes should be set by the user. 
        If ribo == True, lower and upper outliers are removed automatically
    
    """         
    if max_genes == None:
        print("Calculating upper outliers")
        for k in ge.keys():
            adata = adata[~(((adata.obs[groupby]==k) & (adata.obs['n_genes_by_counts']>ge[k])).astype("bool"))]
    else:
        adata = adata[(adata.obs['n_genes_by_counts']<=int(max_genes)).astype("bool"),:].copy()
    print("Now " + str(adata.X.shape[0]) + ' cells after filtering max genes ')
    
    if pct_mt != None:
        adata = adata[adata.obs['pct_counts_mt'] < int(pct_mt)]
    else:
        print("Calculating upper outliers")
        for k in mt.keys():
            adata = adata[~((adata.obs[groupby]==k) & (adata.obs['pct_counts_mt']>mt[k]))]
    print("Now " + str(adata.X.shape[0]) + ' cells after filtering pct mt ')
    
    if total_counts!=None:
        adata = adata[(adata.obs['total_counts']<=int(total_counts)).astype("bool")]
    else:
        print("Calculating upper outliers")
        for k in cn.keys():
            adata = adata[~(((adata.obs[groupby]==k) & (adata.obs['total_counts']>cn[k])).astype("bool"))]
    print("Now " + str(adata.X.shape[0]) + ' cells after filtering total counts ')
    if ribo == True:
        for k in ru.keys():
            adata = adata[~(((adata.obs[groupby]==k) & (adata.obs['pct_counts_ribo']>ru[k])).astype("bool"))]
            adata = adata[~(((adata.obs[groupby]==k) & (adata.obs['pct_counts_ribo']<rl[k])).astype("bool"))]
        print("Now " + str(adata.X.shape[0]) + ' cells after filtering pct ribo ')
    if doublets == True:
        if 'DemuxType_hashing' in adata.obs.columns:
            adata = adata[adata.obs['DemuxType_hashing']=='singlet']
        print("Now " + str(adata.X.shape[0]) + ' cells after keeping singlets by hashing')
        if 'DemuxType_crispr' in adata.obs.columns:
            adata = adata[adata.obs['DemuxType_crispr']=='singlet']
        print("Now " + str(adata.X.shape[0]) + ' cells after keeping singlets by gRNA')
    return adata



def adata_splitter(adata, N_splits):
    cell_batch_size=math.ceil(adata.obs.shape[0]/N_splits)
    batches = []
    for batch_start in range(0, adata.obs.shape[0], cell_batch_size):
        actual_batch_size = min(cell_batch_size, adata.obs.shape[0] - batch_start)
        batch_end = batch_start + actual_batch_size
        partial_df = adata.obs.iloc[batch_start:batch_end, :].index
        partial_adata = adata[(adata.obs.index.isin(partial_df)).astype("bool"),:].copy()
        batches.append(partial_adata)
    return batches
    
    
def centered_embedding_adata(
    adata,
    embeddings_key,
    gene_field,
    guide_field,
    ntc_key,
    pool_ntcs=False
):
    
    """
    Centers an embedding at NTC cells.

    adata: AnnData object
    embeddings_key: key in adata.obsm that contains embedding to be centered, e.g. 'X_pca', 'X_scVI',... Assumed to be a dense
        numpy array. 
    guide_field: field name in adata.obs that contains guide assignment
    gene_field: field name in adata.obs that contains gene-level assignment
    ntc_key: key that corresponds to negative controls in adata.obs[gene_field]
    pool_ntcs: whether to pool NTCs to provide the mean estimate. If False, will first calculate
        means over individual NTC guides and then calculate the overall mean of those.

    Returns the an AnnData object with added centered embeddings_key to adata.obsm.
    """
 
    
    centered_key = embeddings_key+"_centered"
    ntc_idx = adata.obs[gene_field] == ntc_key
    ntc_guides = adata.obs[ntc_idx][guide_field].unique().tolist()
   
    mat = adata.obsm[embeddings_key].copy()
    
    if pool_ntcs:
        mat -= np.mean(mat[ntc_idx], axis=0, keepdims=True)
        
    else:
        means=[]
        for ntc_guide in ntc_guides:
            means.append(np.mean(mat[adata.obs[guide_field].isin(ntc_guides)], axis=0))
        mat -=  np.mean(means, axis=0, keepdims=True)
        
    adata.obsm[centered_key] = mat    
    

def subsample_genes (adata, subsample_dic, subsample_type):
    '''
    A function to randomly subsample cells of specified genes to specific number
    adata: AnnData of our interest
    subsample_dic: a dictionary with genes of interest as keys and subsample size (either as fraction or number of cells)
    subsample_type : a string to show type of subsampling 'n': number of random sample of items  or 'frac': fraction of items to return
    returns a new AnnData with random subset of cells
    '''
    adatas=[]
    genes = [el for el in adata.obs['gene_symbol'].unique()]
    for idx,(key,value) in enumerate(subsample_dic.items()):
        adatas.append(adata[adata.obs["gene_symbol"]==key].copy())
        if subsample_type=='n':
            adatas[idx] = adatas[idx][adatas[idx].obs.sample(n = value).index.copy()].copy()
        elif subsample_type=='frac':
            adatas[idx] = adatas[idx][adatas[idx].obs.sample(frac = value).index.copy()].copy()
        else:
            print("subsample_type should be either n or frac")
        genes.remove(key)
    GENE = adata[adata.obs["gene_symbol"].isin(genes)].copy()
    adatas.append(GENE)
    #bdata = GENE.concatenate(adatas)
    bdata = ad.concat(adatas,join='outer', merge='first')
    bdata.var=bdata.var[["ID", "Symbol","Type"]]
    return bdata


### Sphering
def sphering_transform(cov, reg_param=1e-6, reg_trace=False, gpu=False, rotate=True):
    xp = cp if gpu else np

    # shrink
    if reg_trace:
        cov = (1 - reg_param) * cov + reg_param * xp.trace(cov)/cov.shape[0] * xp.eye(cov.shape[0])
    else:
        cov = cov + reg_param * xp.eye(cov.shape[0])
    s, V = xp.linalg.eigh(cov)

    D = xp.diag(1. / xp.sqrt(s))
    Dinv = xp.diag(xp.sqrt(s))
    W = V @ D
    Winv = Dinv @ V.T
    if rotate:
        W = W @ V.T
        Winv = V @ Winv

    return W, Winv

def empirical_covariance(X, gpu=False, ddof=1.0):
    xp = cp if gpu else np
    # Center data
    n_obs, _ = X.shape
    loc = xp.mean(X, axis=0)
    X = X - loc
    cov = (X.T @ X) / (n_obs - ddof)
    return loc, cov

def oas_covariance(X, gpu=False, ddof=1.0):
    xp = cp if gpu else np
    X = xp.asarray(X)

    # Calculate covariance matrix
    n_obs, n_feat = X.shape
    loc = xp.mean(X, axis=0)
    X = X - loc
    cov = (X.T @ X) / (n_obs - ddof)

    # Calculate sufficient statistics
    tr = xp.trace(cov)
    tr_sq = tr**2
    frob_sq = xp.sum(cov**2)

    # Calculate OAS statistics
    num = (1 - 2/n_feat) * frob_sq + tr_sq
    denom = (n_obs + 1 - 2/n_feat) * (frob_sq - tr_sq/n_feat)
    shrinkage = 1.0 if denom == 0 else min(num/denom, 1.0)
    cov_ret = (1.0 - shrinkage) * cov + shrinkage * tr/n_feat

    return loc, cov_ret

class SpheringTransform(object):
    def __init__(self, controls, reg_param=1e-6, reg_trace=False, rotate=True, gpu=False, oas=False):
        self.gpu = gpu
        
        if self.gpu:
            controls = cp.asarray(controls)
        if oas:
            self.mu, cov = oas_covariance(controls, gpu=gpu)
        else:
            self.mu, cov = empirical_covariance(controls, gpu=gpu)
        self.W, self.Winv = (
            sphering_transform(cov, reg_param=reg_param, reg_trace=reg_trace, gpu=gpu, rotate=rotate)
        )

    def normalize(self, X):
        xp = cp if self.gpu else np
        if self.gpu:
            X = cp.asarray(X)

        ret = (X - self.mu) @ self.W
        return ret

    def recolor(self, X):
        xp = cp if self.gpu else np
        if self.gpu:
            X = cp.asarray(X)

        ret = X @ self.Winv + self.mu
        return ret

def run_sphering_transform(
    adata: ad.AnnData,
    reg_param: float = 1e-3,
    reg_trace: bool = True,
    gpu: bool = False,
    oas: bool = True,
    query: str | None = "gene_symbol == 'NTC'",
    embedding_key: str = "X_pca",
    out_key: str = "X_pca_sphered"
) -> ad.AnnData: 
    if query is not None:
        idx = adata.obs.query(query).index
    else:
        idx = adata.index
    
    control_features = adata[idx,:].obsm[embedding_key]
    if gpu:
        control_features = cp.asarray(control_features)
    
    normalizer = SpheringTransform(
        control_features, reg_param=reg_param, reg_trace=reg_trace, gpu=gpu, oas=oas
    )

    # Apply the normalization transformation
    features = adata.obsm[embedding_key]
    if gpu:
        features = cp.asarray(features)
    features_sphered = normalizer.normalize(adata.obsm[embedding_key])
    if gpu:
        features_sphered = features_sphered.get()

    adata.obsm[out_key] = features_sphered
    
    return adata



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


