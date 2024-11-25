from typing import *
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_fscore_support
from itertools import permutations
from tqdm.autonotebook import tqdm


def parse_geneset(filepath) -> Dict[str, Tuple[str, ...]]:
    with open(filepath, "r") as fh:
        txt = fh.read()
    clusters = {}
    for l in txt.strip().split("\n"):
        n, g = l.strip().split("\t\t")
        clusters[n.strip()] = tuple(i.strip() for i in g.strip().split("\t"))
    return clusters


def adjacency_coords_offdiagonal_triangle(df: pd.DataFrame, lower: bool=False) -> pd.DataFrame:
    """Filter adjacency coordinates to include only rows corresponding to upper/lower offdiagonal triangle
    """
    if lower:
        return df.query("x < y")
    return df.query("x > y")


def adjacency_coords(clusters: Dict[str, Tuple[str, ...]], cm: pd.DataFrame, lower: Optional[bool]=None) -> pd.DataFrame:
    """Generate adjacency coordinates from clusters
    """
    if isinstance(clusters, (list, tuple)):
        assert isinstance(clusters[0], str)
        clusters = {"": clusters}
    out = pd.DataFrame(columns=["Cluster", "Gene (x)", "Gene (y)"])
    index = cm.index.to_list()
    for cluster, genes in clusters.items():
        for i, j in permutations([i for i in genes if i in index], 2):
            out.loc[out.shape[0]] = [cluster, i, j]
    for d in "xy":
        out[d] = [index.index(i) for i in out[f"Gene ({d})"]]
    if lower is not None:
        out = adjacency_coords_offdiagonal_triangle(out, lower)
    return out


def adjacency_matrix(clusters: Union[Dict[str, Any], pd.DataFrame], cm: pd.DataFrame) -> pd.DataFrame:
    """Get dense adjacency matrix from adjacency coordinates
    """
    if isinstance(clusters, dict):
        clusters = adjacency_coords(clusters, cm)
    out = pd.DataFrame(np.zeros_like(cm, dtype=bool), index=cm.index, columns=cm.columns)
    for _, row in clusters.iterrows():
        out.loc[row["Gene (y)"], row["Gene (x)"]] = True
    return out


def offdiagonal_triangle_mask(x: Union[np.ndarray, pd.DataFrame], lower: bool=False, incl_diag: bool=False) -> np.ndarray:
    """Create mask for upper/lower offdiagonal triangle of adjacency matrix
    """
    if isinstance(x, pd.DataFrame):
        x = x.to_numpy()
    y = np.zeros_like(x, dtype=bool)
    if lower:
        for r in range(y.shape[0]):
            y[r, :r] = True
    else:
        for r in range(y.shape[0]):
            y[r, r:] = True
    if not incl_diag:
        y[np.eye(y.shape[0], dtype=bool)] = False
    return y


def offdiagonal_triangle(x: Union[np.ndarray, pd.DataFrame], lower: bool=False, incl_diag: bool=False) -> Union[np.ndarray, pd.DataFrame]:
    """Index into NDArray/DataFrame
    """
    if isinstance(x, pd.DataFrame):
        x = x.to_numpy()
    return x[offdiagonal_triangle_mask(x, lower=lower, incl_diag=incl_diag)]


def threshold_rows_by_percentile(df: Union[np.ndarray, pd.DataFrame], p: float) -> Union[np.ndarray, pd.DataFrame]:
    """Emulate searching within k% most similar genes by row-wise thresholding adjacency matrix
    """
    if isinstance(df, pd.DataFrame):
        x = df.to_numpy().copy()
    else:
        x = df.copy()
    x[np.eye(x.shape[0], dtype=bool)] = -1
    t = np.percentile(x, p, axis=1).reshape(-1, 1)
    y = x >= t
    if isinstance(df, pd.DataFrame):
        y = pd.DataFrame(y, index=df.index, columns=df.columns)
    return y


def precision_recall(cm: pd.DataFrame, adj: pd.DataFrame, percentiles=np.linspace(90, 99, 10).astype(int), relative: bool=True, excl: bool=True, verbose: bool=False):
    assert adj.shape[0] > adj.shape[1]
    adj_mat = adjacency_matrix(adj, cm)

    if verbose:
        adj = adjacency_coords_offdiagonal_triangle(adj)
        assert adj.shape[0] > 0
        adj = adj.copy()
        adj["similarity"] = [cm.loc[row["Gene (y)"], row["Gene (x)"]] for _, row in adj.iterrows()]

    if excl:
        sel = adj_mat.max(0) > 0
        adj_mat = adj_mat.loc[sel, sel]

    triangle = offdiagonal_triangle_mask(adj_mat)
    adj_mat_masked = adj_mat.to_numpy()[triangle]

    results = pd.DataFrame(columns=["percentile", "precision", "recall", "f1", "total_pairs", "total_positives", "total_negatives", "true_positives", "false_positives", "true_negatives", "false_negatives"])
    for percentile in tqdm(percentiles):
        if relative:
            tcm = threshold_rows_by_percentile(cm, percentile)
        else:
            # t = cm.to_numpy().copy()[triangle]
            t = np.percentile(offdiagonal_triangle(cm), percentile)
            # print(t)
            tcm = cm >= t
        if excl:
            tcm = tcm.loc[sel, sel]
        tcm_masked = tcm.to_numpy()[triangle]
        metrics = precision_recall_fscore_support(adj_mat_masked, tcm_masked, average="binary")
        assert metrics[-1] is None

        total_pairs = len(adj_mat_masked)
        total_positives = int(np.sum(adj_mat_masked))
        total_negatives = int(np.sum(adj_mat_masked == False))

        true_positives = int(np.array([tcm_masked, adj_mat_masked]).all(0).sum())
        false_positives = int(np.array([tcm_masked, adj_mat_masked == False]).all(0).sum())

        true_negatives = int(np.array([tcm_masked == False, adj_mat_masked == False]).all(0).sum())
        false_negatives = int(np.array([tcm_masked == False, adj_mat_masked]).all(0).sum())

        results.loc[results.shape[0]] = [percentile, *metrics[:-1], total_pairs, total_positives, total_negatives, true_positives, false_positives, true_negatives, false_negatives]

        if verbose:
            recall = []
            for _, row in adj.iterrows():
                try:
                    r = tcm.loc[row["Gene (y)"], row["Gene (x)"]]
                except KeyError:
                    r = pd.NA
                recall.append(r)
            adj[f"recall_{percentile}"] = recall

    if verbose:
        return results, adj
    return results
