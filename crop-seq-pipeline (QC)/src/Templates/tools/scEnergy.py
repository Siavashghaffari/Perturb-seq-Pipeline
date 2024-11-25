from functools import partial
from typing import Dict

from sklearn.metrics import pairwise_distances
import numpy as np
from dask.distributed import Client
from toolz import valmap


def energy_distance(x: np.ndarray, y: np.ndarray):
    dxx = 0.5 * pairwise_distances(x).mean()
    dyy = 0.5 * pairwise_distances(y).mean()
    dxy = 0.5 * pairwise_distances(x, y).mean()
    return 2 * dxy - dxx - dyy


def energy_distance_test(
    x: np.ndarray, y: np.ndarray, nbootstrap: int = 1000, return_null: bool = False
):
    def e_dist(ix1, ix2):
        return (
            dzz[np.ix_(ix1, ix2)].mean()
            - 0.5 * dzz[np.ix_(ix1, ix1)].mean()
            - 0.5 * dzz[np.ix_(ix2, ix2)].mean()
        )

    z = np.concatenate([x, y])
    dzz = pairwise_distances(z)
    stat = e_dist(*np.split(np.arange(len(z)), [len(x)]))
    null = np.array(
        [
            e_dist(*np.split(np.random.permutation(len(z)), [len(x)]))
            for _ in range(nbootstrap)
        ]
    )
    result = {"pvalue": np.mean(null > stat), "statistic": stat}
    if return_null:
        result["null"] = null
    return result


def run(
    client: Client,
    perturbations: Dict[str, np.ndarray],
    controls: Dict[str, np.ndarray],
    pool: str,
):
    def find_closest_control(x):
        if pool == 'No':
            distances = valmap(partial(energy_distance, x), controls)
            closest = min(distances, key=distances.get)
            return closest, controls[closest]
        else:
            closest = 'pool'
            return closest, controls['pool']

    def worker(perturb_emb):
        control_key, control_emb = find_closest_control(perturb_emb)
        result = energy_distance_test(perturb_emb, control_emb)
        result["control_key"] = control_key
        return result

    futures = {
        key: client.submit(worker, future)
        for key, future in client.scatter(perturbations).items()
    }
    return client.gather(futures)

def get_values_umap(adata, col, value):
    return (adata.obs[col] == value).values