 #basic libraries
import numpy as np
import pandas as pd
import scipy
import os
import sys
from glob import glob
from tqdm import tqdm
from scipy.stats import zscore

#plotting
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns

# single cell
import scanpy as sc
import anndata
import pySingleCellNet as pySCN
#DatasetDB

import pydsdb
#from pydsdb.core import list_experiments
#from pydsdb.core import get_dataset
pd.set_option('display.max_columns', None)

import warnings
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import confusion_matrix

# random_forest replace original sc_makeClassifier function in pySCN code
def add_randCategory(expTrain, genes, groups, nRand=70):
    #This code exists in pySCN, just separating it for clarity.
    randDat = pySCN.randomize(expTrain, num=nRand)
    expT = pd.concat([expTrain, randDat])
    allgenes = expT.columns.values
    missingGenes = np.setdiff1d(np.unique(genes), allgenes)
    ggenes= np.intersect1d(np.unique(genes), allgenes)
    ggroups=np.append(np.array(groups), np.repeat("rand", nRand)).flatten()
    return expT, ggenes, ggroups


def random_forest(expTrain, genes, groups, ntrees=1000, stratify=False):
    #original code had the stratification option. Leave that in, but I've tested with my data and balanced is a good default.
    if not stratify:
        print('Running stratified random forest')
        clf = RandomForestClassifier(n_estimators=ntrees, random_state=100)
    else:
        print('Running balanced random forest')
        clf = RandomForestClassifier(n_estimators=ntrees, class_weight="balanced", random_state=100)    
    #Adding calibration
    clf_isotonic = CalibratedClassifierCV(clf, method='isotonic')
    clf_isotonic.fit(expTrain.loc[:,genes].to_numpy(), groups)
    return clf_isotonic

def preparexclassifier(aTrain, dLevel, nTopGenes = 100, nTopGenePairs = 100,
                   counts_per_cell_after=1e4, scaleMax=10, limitToHVG=False, 
                   normalization = False, include_all_genes = False):
    
    warnings.filterwarnings('ignore')
    stTrain= aTrain.obs
    
    expRaw = aTrain.to_df()
    expRaw = expRaw.loc[stTrain.index.values]

    adNorm = aTrain.copy()
    if normalization == True:
        print("Entered Matrix normalization")
        sc.pp.normalize_per_cell(adNorm, counts_per_cell_after=counts_per_cell_after)
        sc.pp.log1p(adNorm)

        print("HVG")
        if limitToHVG:
            sc.pp.highly_variable_genes(adNorm, min_mean=0.0125, max_mean=4, min_disp=0.5)
            adNorm = adNorm[:, adNorm.var.highly_variable]

        sc.pp.scale(adNorm, max_value=scaleMax)

    expTnorm = adNorm.to_df()
    expTnorm=expTnorm.loc[stTrain.index.values]

    if include_all_genes == False:
        cgenesA, grps, cgenes_list = pySCN.findClassyGenes(expTnorm,stTrain, dLevel = dLevel, topX = nTopGenes)
    else: 
        cgenesA = np.array(aTrain.var.index)
        grps = aTrain.obs[dLevel]
        cgenes_list = dict()
        for g in np.unique(grps):
            cgenes_list[g] = cgenesA
    print("There are ", len(cgenesA), " classification genes\n")
    xpairs= pySCN.ptGetTop(expTnorm.loc[:,cgenesA], grps, cgenes_list, topX=nTopGenePairs, sliceSize=5000)
    print("There are", len(xpairs), "top gene pairs\n")
    pdTrain = pySCN.query_transform(expRaw.loc[:,cgenesA], xpairs)
    print("Finished pair transforming the data\n")
    # Add random category
    expT, ggenes, ggroups = add_randCategory(pdTrain, xpairs, grps, nRand=70)
    return cgenesA, ggenes, ggroups, expT   
    
def scn_train_rf(xpairs, grps, expT, ntrees, stratify=True):
    return random_forest(expT.loc[:, xpairs], genes=xpairs, groups=grps, ntrees=ntrees, stratify=stratify)  
    
def calc_confusion_matrix(y_true, y_pred, y_labels):
    cm = pd.DataFrame(confusion_matrix(y_true, y_pred, labels=y_labels, normalize='true'), columns=y_labels, index=y_labels)*100
    cm.index.name = 'True labels'
    return cm

def plot_cm(cm, **kwargs):
    sns.clustermap(cm, cmap='coolwarm', row_cluster=False, col_cluster=False, linewidth=1, annot=True, fmt='.1f', **kwargs)

def addLabels_cutoff(adata, probability_cutoff):
    undef_id = 'rand'
    # Add max value as a column
    adata.obs['max_val_classifier'] = np.max(adata.to_df(), axis=1)
    adata.obs[f'SCN_class_proba_{probability_cutoff}'] = adata.obs['SCN_class'].astype(str).copy()
    adata.obs.loc[adata.obs['max_val_classifier']<probability_cutoff, f'SCN_class_proba_{probability_cutoff}'] = undef_id
    return adata

def get_margin(adata):
    return adata.to_df().apply(lambda x:x.max()/np.sort(x)[-2] if np.sort(x)[-2]!=0 else np.nan, axis=1)

def add_label_margin(adata, margin_cutoff):
    undef_id = 'rand'
    # Add max value as a column
    adata.obs['margin'] = get_margin(adata)
    adata.obs[f'SCN_class_margin_{margin_cutoff}'] = adata.obs['SCN_class'].astype(str).copy()
    adata.obs.loc[(adata.obs['margin']<margin_cutoff)|(adata.obs['max_val_classifier']<0.5), f'SCN_class_margin_{margin_cutoff}'] = undef_id
    return adata

def matrix_study_label(adata, study_annotation, annotation):
    cnts = adata.obs.groupby([study_annotation, annotation])[study_annotation].count().reset_index(name='cnt')
    cnts = cnts[cnts['cnt']>0]
    cnts['pct'] = cnts.groupby(study_annotation)['cnt'].transform(lambda x:np.round(100*x/x.sum(), 1))
    cnts = cnts.pivot(index=study_annotation, columns=annotation, values='pct').fillna(0)
    return cnts

