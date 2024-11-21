import hotspot

#single cell libraries
import scanpy as sc

#general
import pandas as pd
import sys, os
import numpy as np
from glob import glob
from tqdm import tqdm
import json
import scipy
from scipy.stats import kde
from scipy.stats import zscore
import pickle

#plotting
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns

#
from gprofiler import GProfiler




def pickling(filename, data_structure):
    with open(filename, 'wb') as p:
        pickle.dump(data_structure, p)


def unpickling(pickle_file):
    with open(pickle_file, 'rb') as f:
        return pickle.load(f)


def plot_module_barplot(module, n_top_genes, ax, results_module):
    df = results_module[module].sort_values(by='Z', ascending=False)[:n_top_genes].reset_index()
    sns.barplot(data=df, x='Z', y=0, color='#4863A0', ax=ax)
    ax.set_title('module ' + str(module))
    ax.set_ylabel('')
    ax.set_xlabel('Zscore')
    sns.despine()
    
def plot_module(adata, module, n_top_genes, results_module):
    df = results_module[module].sort_values(by='Z', ascending=False)[:n_top_genes].reset_index()
    x = df['Z']
    y = df[0]
    fig, ax = plt.subplots(1, 7, figsize=(22, 5), gridspec_kw={'width_ratios': [0.8, 2, 2, 2, 2, 2, 2], 'wspace':0.05})
    plot_module_barplot(module, n_top_genes, ax[0], results_module)
    sc.pl.umap(adata, color=str(module), cmap='RdBu_r', ax=ax[1], size=10, frameon=False, show=False, title=f'score module {module}')
    sc.pl.umap(adata, color=y[0], cmap='inferno', ax=ax[2], size=10, frameon=False, show=False)
    sc.pl.umap(adata, color=y[1], cmap='inferno', ax=ax[3], size=10, frameon=False, show=False)
    sc.pl.umap(adata, color=y[2], cmap='inferno', ax=ax[4], size=10, frameon=False, show=False)
    sc.pl.umap(adata, color=y[3], cmap='inferno', ax=ax[5], size=10, frameon=False, show=False)
    sc.pl.umap(adata, color=y[4], cmap='inferno', ax=ax[6], size=10, frameon=False, show=False)
    
def plot_clustermap(df):
    y = df.shape[0]*0.3
    x = y*0.2
    vmax=-np.log10(1E-10)
    cm = sns.clustermap(df, cmap='coolwarm', vmax=vmax, figsize=(x,y), cbar=False, col_cluster=False, square=False, linewidth=1)
    cm.cax.set_visible(False)
    cm.ax_col_dendrogram.set_visible(False)
    
def dominant_col(df):
    return df.columns[np.argmax(df.values, axis=1)]

def get_margin(df):
    return df.apply(lambda x:x.max()/np.sort(x)[-2] if np.sort(x)[-2]!=0 else np.nan, axis=1)

def get_second_best(df):
    return df.apply(lambda x:x.nlargest(2).index[-1], axis=1)

def consolidate_best_modules(x):
    y = x.split(',')
    y = [int(el) for el in y]
    y.sort()
    y = str(y)
    return y.replace('[', '').replace(']', '')