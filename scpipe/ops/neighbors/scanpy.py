from __future__ import annotations
import anndata as ad
import scanpy as sc
import numpy as np
from ...types import Result

class ScanpyNeighbors:
    """
    Build kNN graph from an embedding matrix using Scanpy.
    Returns connectivities/distances + a neighbors metadata dict.
    """
    def __init__(self, n_neighbors: int = 15, metric: str = "euclidean", random_state: int = 0):
        self.cfg = dict(n_neighbors=n_neighbors, metric=metric, random_state=random_state)

    def build(self, X_embed) -> Result:
        tmp = ad.AnnData(X=np.asarray(X_embed))
        sc.pp.neighbors(tmp, n_neighbors=self.cfg["n_neighbors"], metric=self.cfg["metric"], use_rep=None,
                        random_state=self.cfg["random_state"])

        # Normalize metadata to be robust downstream
        neigh = dict(tmp.uns.get("neighbors", {}))
        neigh["connectivities_key"] = "connectivities"
        if "distances" in tmp.obsp:
            neigh["distances_key"] = "distances"
        params = dict(neigh.get("params", {}))
        params.setdefault("method", "umap")
        params.setdefault("use_rep", None)
        params.setdefault("n_pcs", None)
        params.setdefault("n_neighbors", self.cfg["n_neighbors"])
        params.setdefault("metric", self.cfg["metric"])
        neigh["params"] = params

        return Result(
            outputs={
                "obsp.connectivities": tmp.obsp["connectivities"],
                "obsp.distances": tmp.obsp["distances"],
                "uns.neighbors": neigh,
            },
            state=self.cfg,
        )