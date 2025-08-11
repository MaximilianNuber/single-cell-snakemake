from __future__ import annotations
import anndata as ad
import numpy as np
import scanpy as sc
from ...types import Result

class ScanpyUMAP:
    """
    UMAP from precomputed neighbors graphs.
    Expects connectivities (CSR) and optional distances; uses neighbors metadata if provided.
    """
    def __init__(self, n_components: int = 2, min_dist: float = 0.5, spread: float = 1.0, random_state: int = 0):
        self.cfg = dict(n_components=n_components, min_dist=min_dist, spread=spread, random_state=random_state)

    def embed(self, connectivities, distances=None, *, neighbors_uns: dict | None = None) -> Result:
        n = connectivities.shape[0]
        tmp = ad.AnnData(X=np.zeros((n, 1)))
        tmp.obsp["connectivities"] = connectivities
        if distances is not None:
            tmp.obsp["distances"] = distances

        if neighbors_uns is None:
            neighbors_uns = {
                "connectivities_key": "connectivities",
                **({"distances_key": "distances"} if distances is not None else {}),
                "params": {"method": "umap", "use_rep": None, "n_pcs": None, "metric": "precomputed"},
            }
        else:
            neighbors_uns = dict(neighbors_uns)
            neighbors_uns["connectivities_key"] = "connectivities"
            if distances is not None:
                neighbors_uns["distances_key"] = "distances"
            params = dict(neighbors_uns.get("params", {}))
            params.setdefault("method", "umap")
            params.setdefault("use_rep", None)
            params.setdefault("n_pcs", None)
            neighbors_uns["params"] = params

        tmp.uns["neighbors"] = neighbors_uns
        sc.tl.umap(tmp, **self.cfg)
        return Result(outputs={"obsm.X_umap": tmp.obsm["X_umap"]}, state=self.cfg)