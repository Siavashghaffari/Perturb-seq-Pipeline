#
# Copyright (c) 2020, NVIDIA CORPORATION.
# Modified by Siavash Ghaffari to be consistent with the crop-seq-pipeline

#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import rmm
import cupy as cp
import cudf
import cugraph

import time
import dask
from cuml.dask.common.part_utils import _extract_partitions
from cuml.common import with_cupy_rmm

import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import math
import h5py

from cuml.linear_model import LinearRegression
from statsmodels import robust

import warnings

warnings.filterwarnings("ignore", "Expected ")
warnings.simplefilter("ignore")

from cuml.linear_model import LinearRegression
from cuml.preprocessing import StandardScaler
from cuml.decomposition import PCA

from anndata import AnnData
import scvi
from harmony import harmonize

# Helper Scripts
import Templates.tools.scProc as proc


def scale(normalized, max_value=10):
    """
    Scales matrix to unit variance and clips values

    Parameters
    ----------

    normalized : cupy.ndarray or numpy.ndarray of shape (n_cells, n_genes)
                 Matrix to scale
    max_value : int
                After scaling matrix to unit variance,
                values will be clipped to this number
                of std deviations.

    Return
    ------

    normalized : cupy.ndarray of shape (n_cells, n_genes)
        Dense normalized matrix
    """

    scaled = StandardScaler().fit_transform(normalized)

    return cp.clip(scaled, a_min=-max_value, a_max=max_value)


import h5py
from statsmodels import robust


def _regress_out_chunk(X, y):
    """
    Performs a data_chunk.shape[1] number of local linear regressions,
    replacing the data in the original chunk w/ the regressed result.

    Parameters
    ----------

    X : cupy.ndarray of shape (n_cells, 3)
        Matrix of regressors

    y : cupy.ndarray or cupy.sparse.spmatrix of shape (n_cells,)
        containing a single column of the cellxgene matrix

    Returns
    -------

    dense_mat : cupy.ndarray of shape (n_cells,)
        Adjusted column
    """
    if cp.sparse.issparse(y):
        y = y.todense()

    lr = LinearRegression(fit_intercept=False, output_type="cupy")
    lr.fit(X, y, convert_dtype=True)
    return y.reshape(
        y.shape[0],
    ) - lr.predict(
        X
    ).reshape(y.shape[0])


def normalize_total(csr_arr, target_sum):
    """
    Normalizes rows in matrix so they sum to `target_sum`

    Parameters
    ----------

    csr_arr : cupy.sparse.csr_matrix of shape (n_cells, n_genes)
        Matrix to normalize

    target_sum : int
        Each row will be normalized to sum to this value

    Returns
    -------

    csr_arr : cupy.sparse.csr_arr of shape (n_cells, n_genes)
        Normalized sparse matrix
    """

    mul_kernel = cp.RawKernel(
        r"""
    extern "C" __global__
    void mul_kernel(const int *indptr, float *data, 
                    int nrows, int tsum) {
        int row = blockDim.x * blockIdx.x + threadIdx.x;
        
        if(row >= nrows)
            return;
        
        float scale = 0.0;
        int start_idx = indptr[row];
        int stop_idx = indptr[row+1];

        for(int i = start_idx; i < stop_idx; i++)
            scale += data[i];

        if(scale > 0.0) {
            scale = tsum / scale;
            for(int i = start_idx; i < stop_idx; i++)
                data[i] *= scale;
        }
    }
    """,
        "mul_kernel",
    )

    mul_kernel(
        (math.ceil(csr_arr.shape[0] / 32.0),),
        (32,),
        (csr_arr.indptr, csr_arr.data, csr_arr.shape[0], int(target_sum)),
    )

    return csr_arr


def regress_out(normalized, n_counts, percent_mito, batchsize=100, verbose=False):
    """
    Use linear regression to adjust for the effects of unwanted noise
    and variation.

    Parameters
    ----------

    normalized : cupy.sparse.csc_matrix of shape (n_cells, n_genes)
        The matrix to adjust. The adjustment will be performed over
        the columns.

    n_counts : cupy.ndarray of shape (n_cells,)
        Number of genes for each cell

    percent_mito : cupy.ndarray of shape (n_cells,)
        Percentage of genes that each cell needs to adjust for

    batchsize: Union[int,Literal["all"],None] (default: 100)
        Number of genes that should be processed together.
        If `'all'` all genes will be processed together if `normalized.shape[0]` <100000.
        If `None` each gene will be analysed seperatly.
        Will be ignored if cuML version < 22.12

    verbose : bool
        Print debugging information

    Returns
    -------

    outputs : cupy.ndarray
        Adjusted matrix
    """

    regressors = cp.ones((n_counts.shape[0] * 3)).reshape(
        (n_counts.shape[0], 3), order="F"
    )

    regressors[:, 1] = n_counts
    regressors[:, 2] = percent_mito

    outputs = cp.empty(normalized.shape, dtype=normalized.dtype, order="F")

    # cuML gained support for multi-target regression in version 22.12. This
    # removes the need for a Python for loop and speeds up the code
    # significantly.
    cuml_supports_multi_target = LinearRegression._get_tags()["multioutput"]

    if cuml_supports_multi_target and batchsize:
        if batchsize == "all" and normalized.shape[0] < 100000:
            if cp.sparse.issparse(normalized):
                normalized = normalized.todense()
            X = regressors
            # Use SVD algorithm as this is the only algorithm supported in the
            # multi-target regression. In addition, it is more numerically stable
            # than the default 'eig' algorithm.
            lr = LinearRegression(
                fit_intercept=False, output_type="cupy", algorithm="svd"
            )
            lr.fit(X, normalized, convert_dtype=True)
            outputs[:] = normalized - lr.predict(X)
        else:
            if batchsize == "all":
                batchsize = 100
            n_batches = math.ceil(normalized.shape[1] / batchsize)
            for batch in range(n_batches):
                start_idx = batch * batchsize
                stop_idx = min(batch * batchsize + batchsize, normalized.shape[1])
                if cp.sparse.issparse(normalized):
                    arr_batch = normalized[:, start_idx:stop_idx].todense()
                else:
                    arr_batch = normalized[:, start_idx:stop_idx].copy()
                X = regressors
                lr = LinearRegression(
                    fit_intercept=False, output_type="cupy", algorithm="svd"
                )
                lr.fit(X, arr_batch, convert_dtype=True)
                # Instead of "return y - lr.predict(X), we write to outputs to maintain
                # "F" ordering like in the else branch.
                outputs[:, start_idx:stop_idx] = arr_batch - lr.predict(X)
    else:
        if normalized.shape[0] < 100000 and cp.sparse.issparse(normalized):
            normalized = normalized.todense()
        for i in range(normalized.shape[1]):
            if verbose and i % 500 == 0:
                print("Regressed %s out of %s" % (i, normalized.shape[1]))
            X = regressors
            y = normalized[:, i]
            outputs[:, i] = _regress_out_chunk(X, y)

    return outputs


def filter_cells(
    sparse_gpu_array, min_genes, max_genes, rows_per_batch=10000, barcodes=None
):
    """
    Filter cells that have genes greater than a max number of genes or less than
    a minimum number of genes.

    Parameters
    ----------

    sparse_gpu_array : cupy.sparse.csr_matrix of shape (n_cells, n_genes)
        CSR matrix to filter

    min_genes : int
        Lower bound on number of genes to keep

    max_genes : int
        Upper bound on number of genes to keep

    rows_per_batch : int
        Batch size to use for filtering. This can be adjusted for performance
        to trade-off memory use.

    barcodes : series
        cudf series containing cell barcodes.

    Returns
    -------

    filtered : scipy.sparse.csr_matrix of shape (n_cells, n_genes)
        Matrix on host with filtered cells

    barcodes : If barcodes are provided, also returns a series of
        filtered barcodes.
    """

    n_batches = math.ceil(sparse_gpu_array.shape[0] / rows_per_batch)
    filtered_list = []
    barcodes_batch = None
    for batch in range(n_batches):
        batch_size = rows_per_batch
        start_idx = batch * batch_size
        stop_idx = min(batch * batch_size + batch_size, sparse_gpu_array.shape[0])
        arr_batch = sparse_gpu_array[start_idx:stop_idx]
        if barcodes is not None:
            barcodes_batch = barcodes[start_idx:stop_idx]
        filtered_list.append(
            _filter_cells(
                arr_batch,
                min_genes=min_genes,
                max_genes=max_genes,
                barcodes=barcodes_batch,
            )
        )

    if barcodes is None:
        return scipy.sparse.vstack(filtered_list)
    else:
        filtered_data = [x[0] for x in filtered_list]
        filtered_barcodes = [x[1] for x in filtered_list]
        filtered_barcodes = cudf.concat(filtered_barcodes)
        return scipy.sparse.vstack(filtered_data), filtered_barcodes.reset_index(
            drop=True
        )


def _filter_cells(sparse_gpu_array, min_genes, max_genes, barcodes=None):
    degrees = cp.diff(sparse_gpu_array.indptr)
    query = ((min_genes <= degrees) & (degrees <= max_genes)).ravel()
    query = query.get()
    if barcodes is None:
        return sparse_gpu_array.get()[query]
    else:
        return sparse_gpu_array.get()[query], barcodes[query]


def filter_genes(sparse_gpu_array, genes_idx, min_cells=0):
    """
    Filters out genes that contain less than a specified number of cells

    Parameters
    ----------

    sparse_gpu_array : scipy.sparse.csr_matrix of shape (n_cells, n_genes)
        CSR Matrix to filter

    genes_idx : cudf.Series or pandas.Series of size (n_genes,)
        Current index of genes. These must map to the indices in sparse_gpu_array

    min_cells : int
        Genes containing a number of cells below this value will be filtered
    """
    thr = np.asarray(sparse_gpu_array.sum(axis=0) >= min_cells).ravel()
    filtered_genes = cp.sparse.csr_matrix(sparse_gpu_array[:, thr])
    genes_idx = genes_idx[np.where(thr)[0]]

    return filtered_genes, genes_idx.reset_index(drop=True)


def select_groups(labels, groups_order_subset="all"):
    groups_order = labels.cat.categories
    groups_masks = np.zeros(
        (len(labels.cat.categories), len(labels.cat.codes)), dtype=bool
    )
    for iname, name in enumerate(labels.cat.categories):
        # if the name is not found, fallback to index retrieval
        if labels.cat.categories[iname] in labels.cat.codes:
            mask = labels.cat.categories[iname] == labels.cat.codes
        else:
            mask = iname == labels.cat.codes
        groups_masks[iname] = mask.values
    groups_ids = list(range(len(groups_order)))
    if groups_order_subset != "all":
        groups_ids = []
        for name in groups_order_subset:
            groups_ids.append(np.where(name == labels.cat.categories)[0])
        if len(groups_ids) == 0:
            # fallback to index retrieval
            groups_ids = np.where(
                np.in1d(
                    np.arange(len(labels.cat.categories)).astype(str),
                    np.array(groups_order_subset),
                )
            )[0]
        groups_ids = [groups_id.item() for groups_id in groups_ids]
        if len(groups_ids) > 2:
            groups_ids = np.sort(groups_ids)
        groups_masks = groups_masks[groups_ids]
        groups_order_subset = labels.cat.categories[groups_ids].to_numpy()
    else:
        groups_order_subset = groups_order.to_numpy()
    return groups_order_subset, groups_masks


def rank_genes_groups(
    adata,
    groupby,
    groups="all",
    reference="rest",
    n_genes=None,
    **kwds,
):
    """
    Rank genes for characterizing groups.

    Parameters
    ----------

    adata : adata object

    labels : cudf.Series of size (n_cells,)
        Observations groupings to consider

    var_names : cudf.Series of size (n_genes,)
        Names of genes in X

    groups : Iterable[str] (default: 'all')
        Subset of groups, e.g. ['g1', 'g2', 'g3'], to which comparison
        shall be restricted, or 'all' (default), for all groups.

    reference : str (default: 'rest')
        If 'rest', compare each group to the union of the rest of the group.
        If a group identifier, compare with respect to this group.

    n_genes : int (default: 100)
        The number of genes that appear in the returned tables.
    """

    #### Wherever we see "adata.obs[groupby], we should just replace w/ the groups"

    # for clarity, rename variable
    if groups == "all" or groups == None:
        groups_order = "all"
    elif isinstance(groups, (str, int)):
        raise ValueError("Specify a sequence of groups")
    else:
        groups_order = list(groups)
        if isinstance(groups_order[0], int):
            groups_order = [str(n) for n in groups_order]
        if reference != "rest" and reference not in set(groups_order):
            groups_order += [reference]
    labels = pd.Series(adata.obs[groupby]).reset_index(drop="True")
    if reference != "rest" and reference not in set(labels.cat.categories):
        cats = labels.cat.categories.tolist()
        raise ValueError(
            f"reference = {reference} needs to be one of groupby = {cats}."
        )

    groups_order, groups_masks = select_groups(labels, groups_order)

    original_reference = reference

    X = adata.X
    var_names = adata.var_names

    # for clarity, rename variable
    n_genes_user = n_genes
    # make sure indices are not OoB in case there are less genes than n_genes
    if n_genes == None or n_genes_user > X.shape[1]:
        n_genes_user = X.shape[1]
    # in the following, n_genes is simply another name for the total number of genes

    n_groups = groups_masks.shape[0]
    ns = np.zeros(n_groups, dtype=int)
    for imask, mask in enumerate(groups_masks):
        ns[imask] = np.where(mask)[0].size
    if reference != "rest":
        reference = np.where(groups_order == reference)[0][0]
    reference_indices = cp.arange(X.shape[1], dtype=int)

    rankings_gene_scores = []
    rankings_gene_names = []

    # Perform LogReg

    # if reference is not set, then the groups listed will be compared to the rest
    # if reference is set, then the groups listed will be compared only to the other groups listed
    refname = reference
    from cuml.linear_model import LogisticRegression

    reference = groups_order[0]
    if len(groups) == 1:
        raise Exception("Cannot perform logistic regression on a single cluster.")

    grouping_mask = labels.isin(pd.Series(groups_order))
    grouping = labels.loc[grouping_mask]

    X = X[grouping_mask.values, :]
    # Indexing with a series causes issues, possibly segfault

    grouping_logreg = grouping.cat.codes.to_numpy().astype("float32")
    uniques = np.unique(grouping_logreg)
    for idx, cat in enumerate(uniques):
        grouping_logreg[np.where(grouping_logreg == cat)] = idx

    clf = LogisticRegression(**kwds)
    clf.fit(X, grouping_logreg)
    scores_all = cp.array(clf.coef_)
    if len(groups_order) == scores_all.shape[1]:
        scores_all = scores_all.T
    for igroup, group in enumerate(groups_order):
        if len(groups_order) <= 2:  # binary logistic regression
            scores = scores_all[0]
        else:
            scores = scores_all[igroup]

        partition = cp.argpartition(scores, -n_genes_user)[-n_genes_user:]
        partial_indices = cp.argsort(scores[partition])[::-1]
        global_indices = reference_indices[partition][partial_indices]
        rankings_gene_scores.append(scores[global_indices].get())
        rankings_gene_names.append(var_names[global_indices.get()])
        if len(groups_order) <= 2:
            break

    groups_order_save = [str(g) for g in groups_order]
    if len(groups) == 2:
        groups_order_save = [groups_order_save[0]]

    scores = np.rec.fromarrays(
        [n for n in rankings_gene_scores],
        dtype=[(rn, "float32") for rn in groups_order_save],
    )

    names = np.rec.fromarrays(
        [n for n in rankings_gene_names],
        dtype=[(rn, "U50") for rn in groups_order_save],
    )

    return scores, names, original_reference


def leiden(adata, resolution=1.0):
    """
    Performs Leiden Clustering using cuGraph

    Parameters
    ----------

    adata : annData object with 'neighbors' field.

    resolution : float, optional (default: 1)
        A parameter value controlling the coarseness of the clustering.
        Higher values lead to more clusters.

    """
    # Adjacency graph
    adjacency = adata.obsp["connectivities"]
    offsets = cudf.Series(adjacency.indptr)
    indices = cudf.Series(adjacency.indices)
    g = cugraph.Graph()
    if hasattr(g, "add_adj_list"):
        g.add_adj_list(offsets, indices, None)
    else:
        g.from_cudf_adjlist(offsets, indices, None)

    # Cluster
    leiden_parts, _ = cugraph.leiden(g, resolution=resolution)

    # Format output
    clusters = (
        leiden_parts.to_pandas().sort_values("vertex")[["partition"]].to_numpy().ravel()
    )
    clusters = pd.Categorical(clusters)

    return clusters


@with_cupy_rmm
def sq_sum_csr_matrix(client, csr_matrix, axis=0):
    """
    Implements sum operation for dask array when the backend is cupy sparse csr matrix
    """
    client = dask.distributed.default_client()

    def __sq_sum(x):
        x = x.multiply(x)
        return x.sum(axis=axis)

    parts = client.sync(_extract_partitions, csr_matrix)
    futures = [
        client.submit(__sq_sum, part, workers=[w], pure=False) for w, part in parts
    ]
    objs = []
    for i in range(len(futures)):
        obj = dask.array.from_delayed(
            futures[i], shape=futures[i].result().shape, dtype=cp.float32
        )
        objs.append(obj)
    return dask.array.concatenate(objs, axis=axis).compute().sum(axis=axis)


@with_cupy_rmm
def sum_csr_matrix(client, csr_matrix, axis=0):
    """
    Implements sum operation for dask array when the backend is cupy sparse csr matrix
    """
    client = dask.distributed.default_client()

    def __sum(x):
        return x.sum(axis=axis)

    parts = client.sync(_extract_partitions, csr_matrix)
    futures = [client.submit(__sum, part, workers=[w], pure=False) for w, part in parts]
    objs = []
    for i in range(len(futures)):
        obj = dask.array.from_delayed(
            futures[i], shape=futures[i].result().shape, dtype=cp.float32
        )
        objs.append(obj)
    return dask.array.concatenate(objs, axis=axis).compute().sum(axis=axis)


def read_with_filter(
    client,
    sample_file,
    min_genes=200,
    max_genes=6000,
    min_cells=1,
    num_cells=None,
    batch_size=50000,
    partial_post_processor=None,
):
    """
    Reads an h5ad file and applies cell and geans count filter. Dask Array is
    used allow partitioning the input file. This function supports multi-GPUs.
    """

    # Path in h5 file
    _data = "/X/data"
    _index = "/X/indices"
    _indprt = "/X/indptr"
    _genes = "/var/_index"
    _barcodes = "/obs/_index"

    @dask.delayed
    def _read_partition_to_sparse_matrix(
        sample_file,
        total_cols,
        batch_start,
        batch_end,
        min_genes=200,
        max_genes=6000,
        post_processor=None,
    ):
        with h5py.File(sample_file, "r") as h5f:
            indptrs = h5f[_indprt]
            start_ptr = indptrs[batch_start]
            end_ptr = indptrs[batch_end]

            # Read all things data and index
            sub_data = cp.array(h5f[_data][start_ptr:end_ptr])
            sub_indices = cp.array(h5f[_index][start_ptr:end_ptr])

            # recompute the row pointer for the partial dataset
            sub_indptrs = cp.array(indptrs[batch_start : (batch_end + 1)])
            sub_indptrs = sub_indptrs - sub_indptrs[0]
        start = time.time()

        # Reconstruct partial sparse array
        partial_sparse_array = cp.sparse.csr_matrix(
            (sub_data, sub_indices, sub_indptrs),
            shape=(batch_end - batch_start, total_cols),
        )

        # TODO: Add barcode filtering here.
        degrees = cp.diff(partial_sparse_array.indptr)
        query = (min_genes <= degrees) & (degrees <= max_genes)
        partial_sparse_array = partial_sparse_array[query]

        if post_processor is not None:
            partial_sparse_array = post_processor(partial_sparse_array)

        return partial_sparse_array

    with h5py.File(sample_file, "r") as h5f:
        # Compute the number of cells to read
        indptr = h5f[_indprt]
        genes = cudf.Series(h5f[_genes], dtype=cp.dtype("object"))

        total_cols = genes.shape[0]
        max_cells = indptr.shape[0] - 1
        if num_cells is not None:
            max_cells = num_cells

    dls = []
    for batch_start in range(0, max_cells, batch_size):
        actual_batch_size = min(batch_size, max_cells - batch_start)
        dls.append(
            dask.array.from_delayed(
                (_read_partition_to_sparse_matrix)(
                    sample_file,
                    total_cols,
                    batch_start,
                    batch_start + actual_batch_size,
                    min_genes=min_genes,
                    max_genes=max_genes,
                    post_processor=partial_post_processor,
                ),
                dtype=cp.float32,
                shape=(actual_batch_size, total_cols),
            )
        )

    dask_sparse_arr = dask.array.concatenate(dls)
    dask_sparse_arr = dask_sparse_arr.persist()

    # Filter by genes (i.e. cell count per gene)
    gene_wise_cell_cnt = sum_csr_matrix(client, dask_sparse_arr)
    query = gene_wise_cell_cnt > min_cells

    # Filter genes for var
    genes = genes[query]
    genes = genes.reset_index(drop=True)

    query = cp.where(query == True)[0]
    dask_sparse_arr = dask_sparse_arr[:, query.get()].persist()

    return dask_sparse_arr, genes, query


def highly_variable_genes_filter(client, data_mat, genes, n_top_genes=None):
    if n_top_genes is None:
        n_top_genes = genes.shape[0] // 10

    mean = sum_csr_matrix(client, data_mat, axis=0) / data_mat.shape[0]
    mean[mean == 0] = 1e-12

    mean_sq = sq_sum_csr_matrix(client, data_mat, axis=0) / data_mat.shape[0]
    variance = mean_sq - mean**2
    variance *= data_mat.shape[1] / (data_mat.shape[0] - 1)
    dispersion = variance / mean

    df = pd.DataFrame()
    df["genes"] = genes.to_numpy()
    df["means"] = mean.tolist()
    df["dispersions"] = dispersion.tolist()
    df["mean_bin"] = pd.cut(
        df["means"],
        np.r_[-np.inf, np.percentile(df["means"], np.arange(10, 105, 5)), np.inf],
    )

    disp_grouped = df.groupby("mean_bin")["dispersions"]
    disp_median_bin = disp_grouped.median()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        disp_mad_bin = disp_grouped.apply(robust.mad)
        df["dispersions_norm"] = (
            df["dispersions"].values - disp_median_bin[df["mean_bin"].values].values
        ) / disp_mad_bin[df["mean_bin"].values].values

    dispersion_norm = df["dispersions_norm"].values

    dispersion_norm = dispersion_norm[~np.isnan(dispersion_norm)]
    dispersion_norm[::-1].sort()

    if n_top_genes > df.shape[0]:
        n_top_genes = df.shape[0]

    disp_cut_off = dispersion_norm[n_top_genes - 1]
    vaiable_genes = np.nan_to_num(df["dispersions_norm"].values) >= disp_cut_off

    return vaiable_genes


def _cellranger_hvg(mean, mean_sq, genes, n_cells, n_top_genes):
    mean[mean == 0] = 1e-12
    variance = mean_sq - mean**2
    variance *= len(genes) / (n_cells - 1)
    dispersion = variance / mean

    df = pd.DataFrame()
    # Note - can be replaced with cudf once 'cut' is added in 21.08
    df["genes"] = genes.to_numpy()
    df["means"] = mean.tolist()
    df["dispersions"] = dispersion.tolist()
    df["mean_bin"] = pd.cut(
        df["means"],
        np.r_[-np.inf, np.percentile(df["means"], np.arange(10, 105, 5)), np.inf],
    )

    disp_grouped = df.groupby("mean_bin")["dispersions"]
    disp_median_bin = disp_grouped.median()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        disp_mad_bin = disp_grouped.apply(robust.mad)
        df["dispersions_norm"] = (
            df["dispersions"].values - disp_median_bin[df["mean_bin"].values].values
        ) / disp_mad_bin[df["mean_bin"].values].values

    dispersion_norm = df["dispersions_norm"].values

    dispersion_norm = dispersion_norm[~np.isnan(dispersion_norm)]
    dispersion_norm[::-1].sort()

    if n_top_genes is None:
        n_top_genes = genes.shape[0] // 10

    if n_top_genes > df.shape[0]:
        n_top_genes = df.shape[0]

    disp_cut_off = dispersion_norm[n_top_genes - 1]
    variable_genes = np.nan_to_num(df["dispersions_norm"].values) >= disp_cut_off
    return variable_genes


def highly_variable_genes(sparse_gpu_array, genes, n_top_genes=None):
    """
    Identifies highly variable genes using the 'cellranger' method.

    Parameters
    ----------

    sparse_gpu_array : scipy.sparse.csr_matrix of shape (n_cells, n_genes)

    genes : cudf series containing genes

    n_top_genes : number of variable genes
    """

    n_cells = sparse_gpu_array.shape[0]
    mean = sparse_gpu_array.sum(axis=0).flatten() / n_cells
    mean_sq = (
        sparse_gpu_array.multiply(sparse_gpu_array).sum(axis=0).flatten() / n_cells
    )
    variable_genes = _cellranger_hvg(mean, mean_sq, genes, n_cells, n_top_genes)

    return variable_genes


def load_adata(input_file):
    print("Loading in AnnData object.")
    adata = sc.read(input_file, backed="r+")
    return adata

def Batch_maker(
    adata,
    cell_batch_size=100000
):
    n_cells = adata.obs.shape[0] 
    n_genes = adata.var.shape[0]
    n_cells_filtered = 0
    partial_obs = []
    batches = []
    for batch_start in range(0, n_cells, cell_batch_size):
        print("Processed {} cells..".format(batch_start))
        # Get batch indices
        actual_batch_size = min(cell_batch_size, n_cells - batch_start)
        batch_end = batch_start + actual_batch_size
        partial_sparse_array = cp.sparse.csr_matrix(
            #cp.array(((adata.layers['counts'].astype("float"))[batch_start:batch_end]).toarray()), 
            ((adata.layers['counts'].astype("float32"))[batch_start:batch_end]), #Pass sparse matrix directly 
            shape=(batch_end - batch_start, n_genes),
        )
        partial_obs_df = adata.obs.iloc[batch_start:batch_end, :]
        partial_obs.append(partial_obs_df)
        batches.append(partial_sparse_array)

    return batches, partial_obs




def filter_cells(
    adata,
    min_genes=200,
    max_genes=6000,
    max_cells=None,
    cell_batch_size=100000,
):
    print("Filtering cells.")
    if max_cells is not None:
        n_cells = min(max_cells, adata.obs.shape[0])
    else:
        n_cells = adata.obs.shape[0] 

    n_genes = adata.var.shape[0]
    n_cells_filtered = 0
    gene_counts = cp.zeros(shape=(n_genes,), dtype=cp.dtype("float32"))
    partial_obs = []
    batches = []
    for batch_start in range(0, n_cells, cell_batch_size):
        print("Processed {} cells..".format(batch_start))
        # Get batch indices
        actual_batch_size = min(cell_batch_size, n_cells - batch_start)
        batch_end = batch_start + actual_batch_size
        partial_sparse_array = cp.sparse.csr_matrix(
            #cp.array(((adata.layers['counts'].astype("float"))[batch_start:batch_end]).toarray()), 
            ((adata.layers['counts'].astype("float32"))[batch_start:batch_end]), #Pass sparse matrix directly 
            shape=(batch_end - batch_start, n_genes),
        )
        partial_obs_df = adata.obs.iloc[batch_start:batch_end, :]
        # Filter cells in the batch
        degrees = cp.diff(partial_sparse_array.indptr)
        if max_genes != None:
            query = (min_genes <= degrees) & (degrees <= max_genes)
        elif max_genes ==None:
            query = (min_genes <= degrees)
        n_cells_filtered += sum(query)
        partial_sparse_array = partial_sparse_array[query]
        partial_obs_df = partial_obs_df.loc[cp.asnumpy(query).astype(bool), :]
        partial_obs.append(partial_obs_df)
        batches.append(partial_sparse_array)
        # Update gene count # Passing a minlength to be at laest the same shape of gene_count
        gene_counts += cp.bincount(partial_sparse_array.indices, minlength=gene_counts.shape[0]) 
                

    return batches, gene_counts, partial_obs, n_cells_filtered


def filter_genes_and_normalize(
    adata, var_genes_key, batches, gene_counts, min_cells=1, target_sum=1e4
):
    print("Filtering genes and normalizing data.")
    genes = cudf.Series(adata.var[var_genes_key].values, dtype=cp.dtype("object"))
    gene_query = gene_counts >= min_cells
    genes_filtered = genes[gene_query].reset_index(drop=True)
    for i, partial_sparse_array in enumerate(batches):
        # Filter genes
        partial_sparse_array = partial_sparse_array[:, gene_query]
        # Normalize
        partial_sparse_array = normalize_total(
            partial_sparse_array, target_sum=target_sum
        )
        # Log transform
        batches[i] = partial_sparse_array.log1p()

    return batches, genes_filtered

##########
def filter_genes(
    adata, var_genes_key, batches, gene_counts, min_cells=1
):
    print("Filtering genes.")
    genes = cudf.Series(adata.var.index.values, dtype=cp.dtype("object"))
    gene_query = gene_counts >= min_cells
    genes_filtered = genes[gene_query].reset_index(drop=True)
    for i, partial_sparse_array in enumerate(batches):
        # Filter genes
        partial_sparse_array = partial_sparse_array[:, gene_query]
        batches[i] = partial_sparse_array.copy()
        

    return batches, genes_filtered

def normalize(
    adata, batches, target_sum=1e4
):
    print("normalizing data.")
    
    for i, partial_sparse_array in enumerate(batches):
        # Normalize
        partial_sparse_array = normalize_total(
            partial_sparse_array, target_sum=target_sum
        )
        # Log transform
        batches[i] = partial_sparse_array.log1p()

    return batches

#######

def filter_highly_variable_genes(
    batches, genes_filtered, n_cells_filtered, gene_batch_size=2000, n_top_genes=2000
):
    print("Filtering to highly variable genes.")
    # Batch across genes to calculate gene-wise dispersions
    mean = []
    mean_sq = []
    for batch_start in range(0, len(genes_filtered), gene_batch_size):
        # Get batch indices
        print("Processed {} genes..".format(batch_start))
        actual_batch_size = min(gene_batch_size, len(genes_filtered) - batch_start)
        batch_end = batch_start + actual_batch_size

        partial_sparse_array = cp.sparse.vstack(
            [x[:, batch_start:batch_end] for x in batches]
        )

        # Calculate sum per gene
        partial_mean = partial_sparse_array.sum(axis=0) / partial_sparse_array.shape[0]
        mean.append(partial_mean)

        # Calculate sq sum per gene - can batch across genes
        partial_sparse_array = partial_sparse_array.multiply(partial_sparse_array)
        partial_mean_sq = (
            partial_sparse_array.sum(axis=0) / partial_sparse_array.shape[0]
        )
        mean_sq.append(partial_mean_sq)

    mean = cp.hstack(mean).ravel()
    mean_sq = cp.hstack(mean_sq).ravel()

    variable_genes = _cellranger_hvg(
        mean, mean_sq, genes_filtered, n_cells_filtered, n_top_genes
    )

    return variable_genes

#################################################################################
def finish_filtering_adata(batches, partial_obs, genes_filtered):
    sparse_gpu_array = cp.sparse.vstack(
        [partial_sparse_array for partial_sparse_array in batches]
    )
    genes_filtered = genes_filtered.reset_index(drop=True)
    combined_obs = pd.concat(partial_obs)
    print("Completed preprocessing.")

    return sparse_gpu_array, combined_obs, genes_filtered

def finish_normalizing(batches, partial_obs):
    sparse_gpu_array = cp.sparse.vstack(
        [partial_sparse_array for partial_sparse_array in batches]
    )
    combined_obs = pd.concat(partial_obs)
    print("Completed preprocessing.")

    return sparse_gpu_array, combined_obs
####################################################################################\


def finish_filtering(batches, partial_obs, genes_filtered, variable_genes):
    sparse_gpu_array = cp.sparse.vstack(
        [partial_sparse_array[:, variable_genes] for partial_sparse_array in batches]
    )
    genes_filtered = genes_filtered[variable_genes].reset_index(drop=True)
    combined_obs = pd.concat(partial_obs)
    print("Completed preprocessing.")

    return sparse_gpu_array, combined_obs, genes_filtered


def convert_sparse_to_dense(sparse_gpu_array):
    print("Converting sparse GPU array to dense GPU array.")
    dense_gpu_array = sparse_gpu_array.A
    return dense_gpu_array


def compute_pca(dense_gpu_array, n_components, Whiten=False):
    start = time.time()
    pca_out = PCA(n_components=n_components, whiten=Whiten).fit_transform(dense_gpu_array)
    elapsed = time.time() - start
    print("PCA took {} seconds".format(elapsed))
    return pca_out


def compute_neighbors(adata, n_neighbors=15, knn_n_pcs=50,use_rep="X_pca"):
    start = time.time()
    sc.pp.neighbors(
        adata,
        use_rep=use_rep,
        n_neighbors=n_neighbors,
        n_pcs=knn_n_pcs,
        method="rapids",
        random_state=42
    )
    elapsed = time.time() - start
    print("kNN took {} seconds".format(elapsed))


def compute_umap(adata, umap_min_dist=0.3, umap_spread=1.0):
    start = time.time()
    sc.tl.umap(adata, min_dist=umap_min_dist, spread=umap_spread, method="rapids", random_state=11)
    elapsed = time.time() - start
    print("UMAP took {} seconds".format(elapsed))


def compute_louvain(adata):
    start = time.time()
    sc.tl.louvain(adata, flavor="rapids")
    elapsed = time.time() - start
    print("Louvain clustering took {} seconds".format(elapsed))


def compute_diff_exp(adata, label_name, n_genes_group):
    start = time.time()
    cluster_labels = cudf.Series.from_categorical(adata.obs[label_name].cat)
    scores, names, reference = rank_genes_groups(
        adata, groupby="louvain", n_genes=n_genes_group, groups="all", reference="rest"
    )
    elapsed = time.time() - start
    print("Differential expression took {} seconds".format(elapsed))
    return scores, names, reference


class RapidsSingleCellPipeline:
    def __init__(self, input_file, input_format="anndata", var_genes_key="Symbol"):
        if input_format == "anndata":
            self._adata = input_file
            self.var_genes_key = var_genes_key
        else:
            raise AttributeError("Only AnnData input format is supported as of now.")
        rmm.reinitialize(managed_memory=True)
        cp.cuda.set_allocator(rmm.rmm_cupy_allocator)
        self.dense_gpu_array = None
        self.adata = None
        self.adata_ = None

    def filter_data(
        self,
        groupby,
        min_genes=500,
        max_genes=None,
        max_cells=None,
        cell_batch_size=100000,
        min_cells=3,
        gene_batch_size=1200,
        ribo=True, 
        doublets=True,
        pct_mt=None, 
        total_counts=None
    ):
        start = time.time()
        print("starting with " + str(self._adata.X.shape[0]) + ' cells')
        ge, mt, cn, ru, rl = proc.Calculate_outlier_thresholds(self._adata, groupby, ribo, max_genes, pct_mt, total_counts)
        batches, gene_counts, partial_obs, n_cells_filtered = filter_cells(
            self._adata,
            min_genes=min_genes,
            max_genes=max_genes,
            max_cells=max_cells,
            cell_batch_size=cell_batch_size,
        )
        batches, genes_filtered = filter_genes(
            self._adata,
            self.var_genes_key,
            batches,
            gene_counts,
            min_cells=min_cells
        )
        sparse_gpu_array, combined_obs, genes_filtered = finish_filtering_adata(
            batches, partial_obs, genes_filtered
        )
        self.sparse_gpu_array = sparse_gpu_array
        adata_m = AnnData(sparse_gpu_array.get())
        adata_m.layers['counts']=AnnData(sparse_gpu_array.get()).X
        adata_m.var_names = genes_filtered.to_pandas()
        adata_m.obs = combined_obs
        adata_m.var = self._adata.var[self._adata.var.index.isin(genes_filtered.to_arrow().to_pylist())]
        self.genes_filtered  = genes_filtered
        self.n_cells_filtered = n_cells_filtered
        self.adata_ = proc.filter_adata_gpu(adata_m, groupby, min_cells, min_genes, ge, mt, cn, ru, rl, ribo, doublets, max_genes, pct_mt, total_counts)
        elapsed = time.time() - start
        print("Preprocessing took {} seconds".format(elapsed))
        
        
    def run_preprocessing(
        self,
        min_genes=500,
        max_genes=None,
        max_cells=None,
        cell_batch_size=100000,
        min_cells=3,
        gene_batch_size=1200,
        target_sum=1e4,
        n_top_genes=2000,
    ):
        start = time.time()
        self.adata_ = self._adata.copy()
        batches, gene_counts, partial_obs, n_cells_filtered = filter_cells(
            self.adata_,
            min_genes=min_genes,
            max_genes=max_genes,
            max_cells=max_cells,
            cell_batch_size=cell_batch_size,
        )
        batches, genes_filtered = filter_genes_and_normalize(
            self.adata_,
            self.var_genes_key,
            batches,
            gene_counts,
            min_cells=min_cells,
            target_sum=target_sum,
        )
        variable_genes = filter_highly_variable_genes(
            batches,
            genes_filtered,
            n_cells_filtered,
            gene_batch_size=gene_batch_size,
            n_top_genes=n_top_genes,
        )
        sparse_gpu_array, combined_obs, genes_filtered = finish_filtering(
            batches, partial_obs, genes_filtered, variable_genes
        )
        self.adata = AnnData(sparse_gpu_array.get())
        self.adata.var_names = genes_filtered.to_pandas()
        self.adata.obs = combined_obs
        self.dense_gpu_array = convert_sparse_to_dense(sparse_gpu_array)
        del sparse_gpu_array
        elapsed = time.time() - start
        print("Preprocessing took {} seconds".format(elapsed))


    def run_normalize(self,
        cell_batch_size=100000,
        target_sum=1e4):
        start = time.time()
        batches, partial_obs = Batch_maker(
            self.adata_,
            cell_batch_size=cell_batch_size,
        )
        batches = normalize(
            self.adata_,
            batches,
            target_sum=target_sum
        )
        sparse_gpu_array, combined_obs= finish_normalizing(
            batches, partial_obs
        )
        adata_n = AnnData(sparse_gpu_array.get())
        adata_n.obs = combined_obs
        self.adata = adata_n.copy()
        self.adata.layers["counts"] = self.adata.X
        self.batches = batches
        self.partial_obs = partial_obs
        elapsed = time.time() - start
        print("Preprocessing took {} seconds".format(elapsed))
        
        
    def run_highly_variable_genes(
        self,
        gene_batch_size=1200,
        n_top_genes=2000,
    ):
        start = time.time()
        
        variable_genes = filter_highly_variable_genes(
            self.batches,
            self.genes_filtered,
            self.n_cells_filtered,
            gene_batch_size=gene_batch_size,
            n_top_genes=n_top_genes,
        )
        sparse_gpu_array, combined_obs, genes_filtered = finish_filtering(
            self.batches, self.partial_obs, self.genes_filtered, variable_genes
        )
        self.adata = AnnData(sparse_gpu_array.get())
        self.adata.var_names = genes_filtered.to_pandas()
        self.adata.obs = combined_obs
        self.dense_gpu_array = convert_sparse_to_dense(sparse_gpu_array)
        del sparse_gpu_array
        elapsed = time.time() - start
        print("Preprocessing took {} seconds".format(elapsed))

    
    def run_scale(self, max_value=10):
        start = time.time()
        self.dense_gpu_array = scale(self.dense_gpu_array, max_value=10)
        elapsed = time.time() - start
        print("scaling took {} seconds".format(elapsed))
        
        
        
    
    def run_pca(self, n_components):
        pca_out = compute_pca(self.dense_gpu_array, n_components, Whiten=False)
        self.adata.obsm["X_pca"] = pca_out.get()
        
    def run_pca_whitened(self, n_components):
        pca_out = compute_pca(self.dense_gpu_array, n_components, Whiten=True)
        self.adata.obsm["X_pca_whitened"] = pca_out.get()

    def run_knn(self, n_neighbors=15, knn_n_pcs=50, use_rep="X_pca"):
        if "X_pca" in self.adata.obsm:
            compute_neighbors(self.adata, n_neighbors=n_neighbors, knn_n_pcs=knn_n_pcs, use_rep=use_rep)
        else:
            raise RuntimeError("Need to run PCA first.")

    def run_umap(self, umap_min_dist=0.3, umap_spread=1.0):
        if "neighbors" in self.adata.uns:
            compute_umap(
                self.adata, umap_min_dist=umap_min_dist, umap_spread=umap_spread
            )
        else:
            raise RuntimeError("Need to run kNN first.")

    def run_louvain(self):
        if "neighbors" in self.adata.uns:
            compute_louvain(self.adata)
        else:
            raise RuntimeError("Need to run kNN first.")

    def run_leiden(self):
        if "neighbors" in self.adata.uns:
            start = time.time()
            self.adata.obs["leiden"] = leiden(self.adata)
            elapsed = time.time() - start
            print("Leiden clustering took {} seconds".format(elapsed))
        else:
            raise RuntimeError("Need to run kNN first.")

    def run_diff_exp(self, label_name, n_genes_group):
        if label_name in adata.obs.columns:
            scores, names, reference = compute_diff_exp(
                self.adata, label_name=label_name, n_genes_group=n_genes_group
            )
            adata.uns["rank_genes_groups"] = {}
            adata.uns["rank_genes_groups"]["params"] = dict(
                groupby=label_name, method="logreg", reference=reference, use_raw=False
            )
            adata.uns["rank_genes_groups"]["scores"] = scores
            adata.uns["rank_genes_groups"]["names"] = names
        else:
            raise RuntimeError("Specified label name is not found in obs.")
            
    def run_scvi (self, counts_layer = "counts", batch_key = "DemuxAssignment_hashing", n_hidden=128, n_latent=40):
        start = time.time()
        self.adata = self.adata.to_memory()
        self.adata.layers["counts"] = self.adata.X
        scvi.model.SCVI.setup_anndata(self.adata, layer=counts_layer, batch_key=batch_key)
        vae = scvi.model.SCVI(self.adata, n_hidden=128, n_latent=40)
        vae.train(use_gpu=True)
        self.adata.obsm["X_scVI"] = vae.get_latent_representation()
        elapsed = time.time() - start
        print("Scvi Training took {} seconds".format(elapsed))
        
    def run_harmony (self, batch_key = "DemuxAssignment_hashing"):
        start = time.time()
        Z = harmonize(self.adata.obsm['X_pca'], self.adata.obs, batch_key = batch_key)
        self.adata.obsm['X_harmony'] = Z
        elapsed = time.time() - start
        print("Harmony Training took {} seconds".format(elapsed))
    
    def proc (self, n_top_genes, norm=True, scale=False, regress=False, regress_cell_cycle=False, embedding=True, n_components=50, n_neighbors=10, n_pcs=30, batch_key=None, filtered= True, **kwargs):
        
        if filtered==False:
            self.run_preprocessing()
        else:
            if norm == True:
                print('normalizing')
                self.run_normalize()

            print('selecting highly variable genes')
            self.run_highly_variable_genes(n_top_genes=n_top_genes)
        
    
        if regress==True:
            print("regressing")
            sc.pp.regress_out(self.adata, ['total_counts', 'pct_counts_mt'])
        if regress_cell_cycle==True:
            print("regressing cell cycle")
            sc.pp.regress_out(self.adata, ['S_score', 'G2M_score'])
        if scale==True:
            print("scaling")
            self.run_scale(max_value=10)
        if embedding==True:
            print("computing PCA")
            self.run_pca(n_components)
        if batch_key!=None:
            print('batch correcting')
            self.run_scvi (self, counts_layer = "counts", batch_key = batch_key, n_hidden=128, n_latent=40)
            self.run_knn(n_neighbors=30, knn_n_pcs=n_pcs, use_rep='X_scVI')
        else:
            self.run_knn(n_neighbors=n_neighbors, knn_n_pcs=n_pcs)
        print("computing UMAP")
        self.run_umap(umap_min_dist= 0.3, umap_spread=1.0)
        print("computing Leiden clustering")
        self.run_leiden()
        
        
        
    
                                       
                                       
   
                                       
                                       
                                       
        

