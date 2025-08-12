# scpipe/ops/umap_ops.py
from __future__ import annotations
from typing import Dict, Any
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from ..types import Result

def umap_from_latent(
    X_latent: np.ndarray,
    obs_names,
    *,
    min_dist: float = 0.30,
    n_components: int = 2,
    random_state: int = 0,
) -> Result:
    """UMAP on a latent embedding; returns 2D (or n_components-D) coordinates."""
    A = ad.AnnData(X=np.zeros((len(obs_names), 1)))
    A.obs_names = pd.Index(obs_names, name="obs_names")
    A.obsm["X_latent"] = np.asarray(X_latent)
    sc.pp.neighbors(A, use_rep="X_latent")  # if graph already exists you can skip this; kept minimal here
    sc.tl.umap(A, min_dist=min_dist, n_components=n_components, random_state=random_state)
    return Result(
        kind="umap",
        outputs={"obsm.X_umap": A.obsm["X_umap"]},
        state={"schema": "umap/scanpy@v1", "min_dist": min_dist, "n_components": n_components,
               "random_state": random_state},
        metrics={},
    )
