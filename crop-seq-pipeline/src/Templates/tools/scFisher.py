 #basic libraries
import numpy as np
import pandas as pd
import scipy
import os
import sys
from glob import glob
from tqdm import tqdm
from scipy.stats import zscore
import json
#import cellrank as cr
#import decoupler as dc
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy.stats import fisher_exact

#plotting
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns

# single cell
import scanpy as sc
import anndata
from harmony import harmonize
import pySingleCellNet as pySCN

pd.set_option('display.max_columns', None)

def count_features(df, groupby):
    cnts = df.groupby(groupby)[groupby[0]].count().reset_index(name='cnt')
    # filter cells below 3
    cnts = cnts[cnts['cnt']>3]
    cnts['pct'] = cnts.groupby([groupby[0]])['cnt'].transform(lambda x : np.round(100*x/x.sum(), 1))
    return cnts
def add_guideID(df):
    return df.groupby(['geneKO'])['DemuxAssignment_crispr'].transform(lambda x:[el for el in range(len(x))])


def count_guides_by_category(df, category, x, y):
    if 'gene_symbol' in df.columns:
        mapping_dict = df.set_index('DemuxAssignment_crispr')['gene_symbol'].to_dict()
    df =  count_features(df, ['DemuxAssignment_crispr', category]).sort_values(by='DemuxAssignment_crispr', ascending=True)
    df[category] = df[category].astype(str)
    cat_counts = pd.pivot_table(df, index='DemuxAssignment_crispr', columns=category, values='cnt').reset_index()
    # fillna with 0 in category
    cats = df[category].unique()
    #for c in [x] + [y]:
    try:
        cat_counts[x] = cat_counts[x].fillna(0)
    except:
        cat_counts[x] = 0
    try:
        cat_counts[y] = cat_counts[y].fillna(0)
    except:
        cat_counts[y] = 0
    
    try:
        cat_counts['geneKO'] = cat_counts['DemuxAssignment_crispr'].map(mapping_dict)    
    except:
        cat_counts['geneKO'] = cat_counts['DemuxAssignment_crispr'].apply(lambda x:x.split('_')[0])
    cat_counts['guideID'] = add_guideID(cat_counts)
    cat_counts = cat_counts.set_index('DemuxAssignment_crispr')
    return cat_counts


def totalCountsCategory(df, category, x_category, y_category):
    total_category_counts = count_features(df, [category]).set_index(category)
    # Get counts in each category 
    x_category_total = total_category_counts.loc[x_category]['cnt']
    y_category_total = total_category_counts.loc[y_category]['cnt']
    return x_category_total, y_category_total


def dfForFisher(df, category, x_category, y_category):
    x_category_total, y_category_total = totalCountsCategory(df, category, x_category, y_category)
    cat_counts = count_guides_by_category(df, category, x_category, y_category)
    cat_counts['Total_x -' + x_category ] = x_category_total - cat_counts[x_category]
    cat_counts['Total_y -' + y_category ] = y_category_total - cat_counts[y_category]
    return cat_counts


def prep2x2fisher(df, x_category, y_category, guide):
    data = [[df.loc[guide][x_category], df.loc[guide][y_category]], 
            [df.loc[guide]['Total_x -' + x_category], df.loc[guide]['Total_y -' + y_category]]
           ]
    return data


def runFisher(data):
    odd_ratio, p_value = stats.fisher_exact(data, alternative='two-sided')
    return odd_ratio, p_value, -np.log10(p_value)


def runFisherGuides(dfFisher, x_category, y_category):
    df = pd.DataFrame()
    guides = dfFisher.index
    if 'geneKO' in dfFisher.columns:
        mapping_dict = dfFisher['geneKO'].to_dict()
    for guide in guides:
        cols = ['odd_ratio', 'p_value', '-log_p_value']
        data = pd.DataFrame([
            runFisher(
                prep2x2fisher(dfFisher, x_category, y_category, guide)
            )
        ], columns=cols)
        data['guide'] = guide
        try:
            data['geneKO'] = mapping_dict[guide]
        except:
            data['geneKO'] = guide.split('_')[0]
        df = pd.concat([data, df])
    return df

def shuffleLabel(df, xProp):
    part1 = df.sample(frac = xProp)
    part1['labelRandom'] = 'A'
    rest = df.drop(part1.index)
    rest['labelRandom'] = 'B'
    return pd.concat([part1[['DemuxAssignment_crispr','labelRandom']], rest[['DemuxAssignment_crispr','labelRandom']]]) 

def plot_dot(df, x, y, hue, size, ax, **kwargs):
    sns.scatterplot(data=df, x=x, y=y, hue=hue, size=size, palette='coolwarm', ax=ax, **kwargs)
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False, fontsize='xx-small')
    for tick in ax.get_xticklabels():
        tick.set_rotation(90)