# scpipe/ops/graph_ops.py
from __future__ import annotations
from typing import Dict, Any
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from ..types import Result

def neighbors_from_latent_standard(
    X_latent: np.ndarray,
    obs_names,
    *,
    n_neighbors: int = 15,
    method: str = "umap",       # 'umap'/'gauss' etc. accepted by scanpy
    random_state: int = 0,
) -> Result:
    """Build kNN graph from a latent embedding."""
    A = ad.AnnData(X=np.zeros((len(obs_names), 1)))
    A.obs_names = pd.Index(obs_names, name="obs_names")
    A.obsm["X_latent"] = np.asarray(X_latent)
    sc.pp.neighbors(A, use_rep="X_latent", n_neighbors=n_neighbors, method=method, random_state=random_state)
    return Result(
        kind="neighbors",
        outputs={"obsp.connectivities": A.obsp["connectivities"], "obsp.distances": A.obsp.get("distances")},
        state={"schema": "neighbors/standard@v1", "n_neighbors": n_neighbors, "method": method,
               "random_state": random_state},
        metrics={},
    )

def neighbors_from_latent_bbknn(
    X_latent: np.ndarray,
    obs,
    *,
    batch_key: str = "sample",
    n_neighbors: int = 15,
) -> Result:
    """Batch-balanced kNN using BBKNN over the latent embedding."""
    import scanpy.external as sce
    A = ad.AnnData(X=np.zeros((len(obs), 1)))
    A.obs = pd.DataFrame(obs).copy()
    A.obsm["X_latent"] = np.asarray(X_latent)
    # heuristic: neighbors per batch
    n_batches = int(pd.Series(A.obs[batch_key]).nunique()) or 1
    within = max(3, n_neighbors // n_batches)
    # bbknn supports use_rep; fall back to default if older version
    try:
        sce.pp.bbknn(A, batch_key=batch_key, neighbors_within_batch=within, use_rep="X_latent")
    except TypeError:
        sce.pp.bbknn(A, batch_key=batch_key, neighbors_within_batch=within)
    return Result(
        kind="neighbors",
        outputs={"obsp.connectivities": A.obsp["connectivities"], "obsp.distances": A.obsp.get("distances")},
        state={"schema": "neighbors/bbknn@v1", "batch_key": batch_key, "n_neighbors": n_neighbors,
               "neighbors_within_batch": within},
        metrics={},
    )
