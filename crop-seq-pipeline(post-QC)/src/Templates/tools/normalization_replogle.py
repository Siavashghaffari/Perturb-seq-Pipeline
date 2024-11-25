# Perturbseq library for loading and manipulating single-cell experiments
# Copyright (C) 2019  Thomas Norman

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# Modified by Sandra Melo Carlos Jan 2024

import pandas as pd
import numpy as np
from sklearn import preprocessing as pre
from pandas.api.types import is_numeric_dtype
from time import time
import gc
from tqdm import tqdm_notebook
from scipy.sparse import issparse
from scipy.sparse import csr_matrix 
import anndata as ad
import cupy as cp


def check_matrix_and_densify(input_obj):
    # Check if input_obj is a sparse array
    if issparse(input_obj):
        print("     Densifying matrix...")
        matrix = pd.DataFrame(input_obj.todense())
        return matrix
    else:
        print("This is not a sparse array.")
        
        


def normalize_matrix_to_control(matrix, control_matrix, scale_by_total=True, median_umi_count=None):
    """ Normalize expression distribution relative to a population of control (unperturbed) cells.
    The normalization proceeds by first normalizing the UMI counts within each cell to the median
    UMI count within the population. The control population is normalized to the same amount. 
    The expression within the population is then Z-normalized with respect to the control 
    distribution: i.e., for each gene the control mean is subtracted, and then these values are
    divided by the control standard deviation.
        
    Args:
        matrix: gene expression matrix to normalize (output from cellranger)
        control_matrix: gene expression matrix of control population
        scale_by_total: Rescale UMI counts within each cell to population median (default: True)
        median_umi_count: If provided, set the median to normalize to. This is useful for example
            when normalizing multiple independent lanes/gemgroups.
    
    Returns:
        DataFrame of normalized expression data
        
    Example:
        >>>control_cells = pop.cells[pop.cells['perturbation']=='DMSO_control'].index.values
        >>>pop.normalized_matrix = normalize_expression_to_control(pop.matrix,
                                                                   pop.matrix.loc[control_cells].copy())
    """
    
    matrix = check_matrix_and_densify(matrix)
    control_matrix = check_matrix_and_densify(control_matrix)
    
    if (scale_by_total):
        print("     Determining scale factors...")
        reads_per_bc = matrix.sum(axis=1)
        
        if median_umi_count is None:
            median_reads_per_bc = np.median(reads_per_bc)
        else:
            median_reads_per_bc = median_umi_count
        
        scaling_factors = median_reads_per_bc / reads_per_bc
        
        print("     Normalizing matrix to median")
        m = matrix.astype(np.float32)
        # Normalize expression within each cell by median total count
        m = m.mul(scaling_factors, axis=0)
        if np.mean(median_reads_per_bc) < 5000:
            print("Scaling with a small number of reads. Are you sure this is what you want?")
            
        control_reads_per_bc = control_matrix.sum(axis=1)
        
        print("     Normalizing control matrix to median")
        control_scaling_factors = median_reads_per_bc / control_reads_per_bc
        c_m = control_matrix.astype(np.float32)
        c_m = c_m.mul(control_scaling_factors, axis=0)
    else:
        m = matrix.astype(np.float32)
        c_m = matrix.astype(np.float32)

    control_mean = c_m.mean()
    control_std = c_m.std()
    
    print("     Scaling matrix to control")
    # Now center and rescale the expression of each gene to average 0 and std 1
    m_out = (m - control_mean).div(control_std, axis=1)
    
    print("     Done.")
    return pd.DataFrame(m_out, columns=m.columns, index=m.index)

def normalize_matrix_to_control_adata(adata, control_key="gene_symbol == 'NTC'", scale_by_total=True, median_umi_count=None):
    idx = adata.obs.query(control_key).index
    control_matrix = adata[idx,:].layers['counts']
    matrix = adata.layers['counts']
    normalized = normalize_matrix_to_control(matrix, control_matrix, scale_by_total=scale_by_total, median_umi_count=median_umi_count)
    normalized.index = adata.obs.index
    normalized.columns = adata.var.index
    return normalized


def proc_gem(small_adata, median_umi_count, control_key="gene_symbol == 'NTC'"):
    gem_group_matrices = normalize_matrix_to_control_adata(
                small_adata, 
                control_key=control_key,
                scale_by_total=True, 
                median_umi_count=median_umi_count)
    gem_group_matrices= csr_matrix(gem_group_matrices.astype("float32"))
    gem_group_matrices = ad.AnnData(X=gem_group_matrices)
    gem_group_matrices.obs = small_adata.obs.copy()
    gem_group_matrices.var = small_adata.var.copy()
    return gem_group_matrices