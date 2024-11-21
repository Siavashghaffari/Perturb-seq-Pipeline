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

#DatasetDB
import pydsdb
from pydsdb import get_datasets

# Helper Scripts
import Templates.tools.scProc as proc

def adata_cleaner (adata,adatas, experiment, fix_barcodes, alt_experiments):
    """
    This function cleans the anndata object and make it ready for the analysis
    """
    # create a tmp anndata object
    adata_tmp = adatas[experiment]
    if len(alt_experiments)==2:
        search_result =[col for col in adata.obs.columns if experiment in col]
        df = adata.obs[search_result]
        df.columns = ['demux_type', 'assignment']
        df2 = adata.obs[['Sample', 'Barcode']]
        DF =pd.concat([df2, df], axis=1)
    else:
        DF=adata.obs.copy()
    # transfer DemuxEM info
    adata_tmp.obs = DF.copy()
    # Convert back the tmp anndata object to the anndata 
    adata = adata_tmp.copy()
    # to be consistent with our tradition store layers counts into adata.X
    adata.X = adata.layers['counts']
    
    if fix_barcodes == True:
        valid = proc.read_demuxEM_map(valid_assignments)
        adata.obs['assignment'] = proc.fix_assignment(adata.obs)
        adata.obs['assignment'] = adata.obs['assignment'].map(valid)
        for k,v in valid.items():
            adata.obs.loc[adata.obs['assignment']==v, 'demux_type'] = 'valid_doublet'

    
    return adata