# scpipe/ops/cluster_ops.py
from __future__ import annotations
from typing import Optional, Dict, Any
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix
from ..types import Result

def leiden_from_connectivities(
    connectivities: csr_matrix,
    obs_names,
    *,
    resolution: float = 1.0,
    key: str = "leiden",
) -> Result:
    """Leiden clustering from a kNN graph."""
    A = ad.AnnData(X=np.zeros((len(obs_names), 1)))
    A.obs_names = pd.Index(obs_names, name="obs_names")
    A.obsp["connectivities"] = connectivities
    sc.tl.leiden(A, resolution=resolution, key_added=key)
    labels = A.obs[key].astype(str).to_numpy()
    return Result(
        kind="clusters",
        outputs={"obs.leiden": labels},
        state={"schema": "clusters/leiden@v1", "resolution": resolution, "key": key},
        metrics={},
    )
