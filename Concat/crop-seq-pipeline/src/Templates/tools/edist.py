import pandas as pd
import numpy as np
import scanpy as sc

from tqdm import tqdm
from sklearn.metrics import pairwise_distances
from warnings import warn
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed

def pairwise_pca_distances(adata, obs_key, obsm_key='X_pca', dist='sqeuclidean',
                           sample_correct=True, verbose=True):
    """Average of pairwise PCA distances between cells of each group in obs_key.
    For each pair of groups defined in adata.obs[obs_key] (e.g. perturbations)
    computes all pairwise distances between cells in adata.obsm[obsm_key] (e.g. PCA space)
    and averages them per group-pair. This results in a distance matrix between all groups.
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    obs_key: `str` in adata.obs.keys()
        Key in adata.obs specifying the groups to consider.
    obsm_key: `str` in adata.obsm (default: `adata.obsm['X_pca']`)
        Key for embedding coordinates to use.
    dist: `str` for any distance in scipy.spatial.distance (default: `sqeuclidean`)
        Distance metric to use in embedding space.
    sample_correct: `bool` (default: `True`)
        Whether make the estimator for sigma more unbiased (dividing by N-1 instead of N, similar to sample and population variance).
    verbose: `bool` (default: `True`)
        Whether to show a progress bar iterating over all groups.
    Returns
    -------
    pwd: pandas.DataFrame
        DataFrame with average of pairwise PCA distances between all groups.
    """

    if obsm_key=='X_pca' and 'X_pca' not in adata.obsm.keys():
        warn('PCA embedding not found, computing...')
        sc.pp.pca(adata)
    
    X = adata.obsm[obsm_key].copy()
    y = adata.obs[obs_key].astype(str).copy()
    groups = pd.unique(y)
    df = pd.DataFrame(index=groups, columns=groups, dtype=float)
    fct = tqdm if verbose else lambda x: x
    for i, p1 in enumerate(fct(groups)):
        x1 = X[y==p1].copy()
        N = len(x1)
        for p2 in groups[i:]:
            x2 = X[y==p2].copy()
            pwd = pairwise_distances(x1, x2, metric=dist)
            M = len(x2)-1 if (p1==p2) & sample_correct else len(x2)
            factor = N * M
            mean_pwd = np.sum(pwd) / factor
            df.loc[p1, p2] = mean_pwd
            df.loc[p2, p1] = mean_pwd
    df.index.name = obs_key
    df.columns.name = obs_key
    df.name = 'pairwise PCA distances'
    return df

def edist(adata, obs_key='perturbation', obsm_key='X_pca', pwd=None, 
          dist='sqeuclidean', sample_correct=0, verbose=True):
    """Computes the edistance to control. Accepts precomputed pwd.
    Computes the pairwise E-distances between all groups of cells defined in
    adata.obs[obs_key] (e.g. perturbations). Distances are computed in embedding
    space given by adata.obsm[obsm_key] (e.g. PCA space).
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    obs_key: `str` in adata.obs.keys() (default: `perturbation`)
        Key in adata.obs specifying the groups to consider.
    obsm_key: `str` in adata.obsm (default: `adata.obsm['X_pca']`)
        Key for embedding coordinates to use.
    dist: `str` for any distance in scipy.spatial.distance (default: `sqeuclidean`)
        Distance metric to use in embedding space.
    sample_correct: `bool` (default: `True`)
        Whether make the estimator for sigma more unbiased (dividing by N-1 instead of N, similar to sample and population variance).
    verbose: `bool` (default: `True`)
        Whether to show a progress bar iterating over all groups.
    Returns
    -------
    estats: pandas.DataFrame
        DataFrame with pairwise E-distances between all groups.
    """
    pwd = pairwise_pca_distances(adata, obs_key=obs_key, obsm_key=obsm_key, 
                                 dist=dist, sample_correct=sample_correct, 
                                 verbose=verbose) if pwd is None else pwd
    # derive basic statistics
    sigmas = np.diag(pwd)
    deltas = pwd
    estats = 2 * deltas - sigmas - sigmas[:, np.newaxis]
    return estats

def onesided_pca_distances(adata, obs_key, control, obsm_key='X_pca', 
                           dist='sqeuclidean', sample_correct=True, 
                           verbose=True):
    """Average of pairwise PCA distances between cells of each group in obs_key with control group.
    For each group defined in adata.obs[obs_key] (e.g. perturbations)
    computes all pairwise distances between cells in adata.obsm[obsm_key] (e.g. PCA space)
    and averages them per group-control-pair. This results in a distance vector with a value for each group.
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    obs_key: `str` in adata.obs.keys()
        Key in adata.obs specifying the groups to consider.
    control: `str` of a category in adata.obs[obs_key]
        Group in obs_key for control cells.
    obsm_key: `str` in adata.obsm (default: `adata.obsm['X_pca']`)
        Key for embedding coordinates to use.
    dist: `str` for any distance in scipy.spatial.distance (default: `sqeuclidean`)
        Distance metric to use in embedding space.
    sample_correct: `bool` (default: `True`)
        Whether make the estimator for sigma more unbiased (dividing by N-1 instead of N, similar to sample and population variance).
    verbose: `bool` (default: `True`)
        Whether to show a progress bar iterating over all groups.
    Returns
    -------
    pwd: pandas.DataFrame
        DataFrame with average PCA distances to control for all groups.
    """

    if obsm_key=='X_pca' and 'X_pca' not in adata.obsm.keys():
        warn('PCA embedding not found, computing...')
        sc.pp.pca(adata)

    groups = pd.unique(adata.obs[obs_key])
    assert control in groups, f'No cells of control group "{control}" were not found in groups defined by "{obs_key}".'
    df = pd.DataFrame(index=groups, columns=['distance'], dtype=float)
    fct = tqdm if verbose else lambda x: x
    
    x1 = adata[adata.obs[obs_key]==control].obsm[obsm_key].copy()
    N = len(x1)
    for p in fct(groups):
        x2 = adata[adata.obs[obs_key]==p].obsm[obsm_key].copy()
        pwd = pairwise_distances(x1, x2, metric=dist)
        M = len(x2)-1 if (p==control) & sample_correct else len(x2)
        factor = N * M  # Thanks to Garrett Wong for finding this bug
        mean_pwd = np.sum(pwd) / factor
        df.loc[p] = mean_pwd
    df.index.name = obs_key
    df.name = f'PCA distances to {control}'
    return df

def self_pca_distances(adata, obs_key, obsm_key='X_pca', dist='sqeuclidean', 
                       sample_correct=True, verbose=True):
    """Average of pairwise PCA distances between cells within each group in obs_key.
    For each group defined in adata.obs[obs_key] (e.g. perturbations)
    computes all pairwise distances between cells within in the space given by adata.obsm[obsm_key] (e.g. PCA space)
    and averages them per group. This results in a distance vector with a value for each group.
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    obs_key: `str` in adata.obs.keys()
        Key in adata.obs specifying the groups to consider.
    obsm_key: `str` in adata.obsm (default: `adata.obsm['X_pca']`)
        Key for embedding coordinates to use.
    dist: `str` for any distance in scipy.spatial.distance (default: `sqeuclidean`)
        Distance metric to use in embedding space.
    sample_correct: `bool` (default: `True`)
        Whether make the estimator for sigma more unbiased (dividing by N-1 instead of N, similar to sample and population variance).
    verbose: `bool` (default: `True`)
        Whether to show a progress bar iterating over all groups.
    Returns
    -------
    pwd: pandas.DataFrame
        DataFrame with average PCA distances to control for all groups.
    """

    if obsm_key=='X_pca' and 'X_pca' not in adata.obsm.keys():
        warn('PCA embedding not found, computing...')
        sc.pp.pca(adata)

    groups = pd.unique(adata.obs[obs_key])
    df = pd.DataFrame(index=groups, columns=['distance'], dtype=float)
    fct = tqdm if verbose else lambda x: x
    
    for p in fct(groups):
        x = adata[adata.obs[obs_key]==p].obsm[obsm_key].copy()
        pwd = pairwise_distances(x, x, metric=dist)
        N = len(x)
        factor = N*(N-1) if sample_correct else N**2
        mean_pwd = np.sum(pwd) / factor
        df.loc[p] = mean_pwd
    df.index.name = obs_key
    df.name = 'PCA distances within groups'
    return df

def edist_to_control(adata, obs_key='perturbation', control='control', 
                     obsm_key='X_pca', dist='sqeuclidean',  
                     sample_correct=True, verbose=True):
    """Computes the edistance to control.
    Computes the all E-distances between all groups of cells defined in
    adata.obs[obs_key] (e.g. perturbations) and control cells. Distances are computed in embedding
    space given by adata.obsm[obsm_key] (e.g. PCA space).
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    obs_key: `str` in adata.obs.keys() (default: `perturbation`)
        Key in adata.obs specifying the groups to consider.
    obsm_key: `str` in adata.obsm (default: `adata.obsm['X_pca']`)
        Key for embedding coordinates to use.
    dist: `str` for any distance in scipy.spatial.distance (default: `sqeuclidean`)
        Distance metric to use in embedding space.
    sample_correct: `bool` (default: `True`)
        Whether make the estimator for sigma more unbiased (dividing by N-1 instead of N, similar to sample and population variance).
    verbose: `bool` (default: `True`)
        Whether to show a progress bar iterating over all groups.
    Returns
    -------
    ed_to_c: pandas.DataFrame
        DataFrame with E-distances between all groups and control group.
    """
    deltas_to_c = onesided_pca_distances(adata, obs_key=obs_key, control=control, 
                                         obsm_key=obsm_key, dist=dist, 
                                         sample_correct=sample_correct, 
                                         verbose=verbose)
    sigmas = self_pca_distances(adata, obs_key, obsm_key=obsm_key, dist=dist, 
                                sample_correct=sample_correct, verbose=False)
    # derive basic statistics
    ed_to_c = 2 * deltas_to_c - sigmas - sigmas.loc[control]
    return ed_to_c









# TODO make etest allow for multiple controls (accept list of controls)

def etest(adata, obs_key='perturbation', obsm_key='X_pca', dist='sqeuclidean',
          control='control', alpha=0.05, runs=1000, sample_correct=True, n_jobs=1,
          correction_method='holm-sidak', verbose=True):
    """Performs Monte Carlo permutation test with E-distance as test statistic.
    Tests for each group of cells defined in adata.obs[obs_key] if it is significantly
    different from control based on the E-distance in adata.obsm[obsm_key] space.
    Does multiple-testing correction using per default with Holm-Sidak.
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    obs_key: `str` in adata.obs.keys() (default: `perturbation`)
        Key in adata.obs specifying the groups to consider.
    obsm_key: `str` in adata.obsm (default: `adata.obsm['X_pca']`)
        Key for embedding coordinates to use.
    dist: `str` for any distance in scipy.spatial.distance (default: `sqeuclidean`)
        Distance metric to use in embedding space.
    control: `str` (default: `'control'`)
        Defines the control group in adata.obs[obs_key] to test against.
    alpha: `float` between `0` and `1` (default: `0.05`)
        significance cut-off for the test to annotate significance.
    runs: `int` (default: `100`)
        Number of iterations for the permutation test. This is basically the resolution of the E-test p-value.
        E.g. if you choose two iterations, then the p-value can only have 3 values `(0, .5, 1)`. Lower numbers will be much faster.
        We do not recommend going lower than `100` and suggest between `100` and `10000` iterations.
    sample_correct: `bool` (default: `True`)
        Whether make the estimator for sigma more unbiased (dividing by N-1 instead of N, similar to sample and population variance).
    n_jobs: `int` (default: `1`)
        Number of jobs to use for parallelization. If `n_jobs=1`, no parallelization is used.
    correction_method: `None` or any valid method for statsmodels.stats.multitest.multipletests (default: `'holm-sidak'`)
        Method used for multiple-testing correction, since we are testing each group in `adata.obs[obs_key]`.
    verbose: `bool` (default: `True`)
        Whether to show a progress bar iterating over all groups.
    Returns
    -------
    tab: pandas.DataFrame
        E-test results for each group in adata.obs[obs_key] with columns
        - edist: E-distance to control
        - pvalue: E-test p-value if group is different from control
        - significant: If p-value < alpha
        - pvalue_adj: Multiple-testing corrected E-test p-value
        - significant_adj: If p-value_adj < alpha
    """

    groups = pd.unique(adata.obs[obs_key])
    
    # Precompute pairwise distances selectively once
    # (we need pairwise distances within each group and between each group and control)
    # Note: this could be improved further, since we compute distances within control multiple times here. Speedup likely minimal though.
    pwds = {}
    for group in groups:
        x = adata[adata.obs[obs_key].isin([group, control])].obsm[obsm_key].copy()
        pwd = pairwise_distances(x, x, metric=dist)
        pwds[group] = pwd

    # Approximate sampling from null distribution (equal distributions)
    fct = tqdm if verbose else lambda x: x  # progress bar y/n
    M = np.sum(adata.obs[obs_key]==control)  # number of cells in control group
    def one_step():
        # per perturbation, shuffle with control and compute e-distance
        df = pd.DataFrame(index=groups, columns=['edist'], dtype=float)
        for group in groups:
            if group==control:
                # Nothing to test here
                df.loc[group] = [0]
                continue
            N = np.sum(adata.obs[obs_key]==group)
            # shuffle the labels
            labels = adata.obs[obs_key].values[adata.obs[obs_key].isin([group, control])]
            shuffled_labels = np.random.permutation(labels)

            # use precomputed pairwise distances
            sc_pwd = pwds[group]  # precomputed pairwise distances between single cells
            idx = shuffled_labels==group

            # Note that this is wrong: sc_pwd[idx, ~idx] but this is correct: sc_pwd[idx, :][:, ~idx]
            # The first produces a vector, the second a matrix (we need the matrix)
            factor = N / (N-1) if sample_correct else 1
            factor_c = M / (M-1) if sample_correct else 1
            delta = np.sum(sc_pwd[idx, :][:, ~idx]) / (N * M)
            sigma = np.sum(sc_pwd[idx, :][:, idx]) / (N * N) * factor
            sigma_c = np.sum(sc_pwd[~idx, :][:, ~idx]) / (M * M) * factor_c

            edistance = 2 * delta - sigma - sigma_c

            df.loc[group] = edistance
        return df.sort_index()
    res = Parallel(n_jobs=n_jobs)(delayed(one_step)() for i in fct(range(runs)))
    
    # the following is faster than the above and produces the same result
    original = []
    for group in groups:
        if group==control:
            # Nothing to test here
            original.append(0)
            continue
        N = np.sum(adata.obs[obs_key]==group)
        # do *not* shuffle the labels
        labels = adata.obs[obs_key].values[adata.obs[obs_key].isin([group, control])]
        
        # use precomputed pairwise distances
        sc_pwd = pwds[group]  # precomputed pairwise distances between single cells
        idx = labels==group
        
        # Note that this is wrong: sc_pwd[idx, ~idx] but this is correct: sc_pwd[idx, :][:, ~idx]
        # The first produces a vector, the second a matrix (we need the matrix)
        factor = N / (N-1) if sample_correct else 1
        factor_c = M / (M-1) if sample_correct else 1
        delta = np.sum(sc_pwd[idx, :][:, ~idx]) / (N * M)
        sigma = np.sum(sc_pwd[idx, :][:, idx]) / (N * N) * factor
        sigma_c = np.sum(sc_pwd[~idx, :][:, ~idx]) / (M * M) * factor_c
        
        edistance = 2 * delta - sigma - sigma_c
        original.append(edistance)
    df = pd.DataFrame(original, index=groups, columns=['edist'])
    df = df.sort_index()

    # Evaluate test (hypothesis vs null hypothesis)
    # count times shuffling resulted in larger or equal e-distance
    results = np.array(pd.concat([r['edist'] - df['edist'] for r in res], axis=1) >= 0, dtype=int)
    n_failures = pd.Series(np.clip(np.sum(results, axis=1), 1, np.inf), index=df.index)
    pvalues = n_failures / runs  # chance that our results was obtained by chance

    # Apply multiple testing correction
    significant_adj, pvalue_adj, _, _ = multipletests(pvalues.values, alpha=alpha, method=correction_method)

    # Aggregate results
    tab = pd.DataFrame({'edist': df['edist'], 'pvalue': pvalues, 
                        'significant': pvalues < alpha, 'pvalue_adj': pvalue_adj, 
                        'significant_adj': significant_adj}, index=df.index)
    return tab



def equal_subsampling(adata, obs_key, N_min=None):
    """Subsample cells while retaining same class sizes.
    Classes are given by obs_key pointing to categorical in adata.obs.
    If N_min is given, downsamples to at least this number instead of the number
    of cells in the smallest class and throws out classes with less than N_min cells.
    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    obs_key: `str` in adata.obs.keys() (default: `perturbation`)
        Key in adata.obs specifying the groups to consider.
    N_min: `int` or `None` (default: `None`)
        If N_min is given, downsamples to at least this number instead of the number
        of cells in the smallest class and throws out classes with less than N_min cells.
    Returns
    -------
    subdata: :class:`~anndata.AnnData`
        Subsampled version of the original annotated data matrix.
    """

    counts = adata.obs[obs_key].value_counts()
    if N_min is not None:
        groups = counts.index[counts>=N_min]  # ignore groups with less than N_min cells to begin with
    else:
        groups=counts.index
    # We select downsampling target counts by min-max, i.e.
    # the largest N such that every group has at least N cells before downsampling.
    N = np.min(counts)
    N = N if N_min==None else np.max([N_min, N])
    # subsample indices per group
    indices = [np.random.choice(adata.obs_names[adata.obs[obs_key]==group], size=N, replace=False) for group in groups]
    selection = np.hstack(np.array(indices))
    return adata[selection].copy()