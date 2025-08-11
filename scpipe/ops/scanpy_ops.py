from __future__ import annotations
# lib/ops/scanpy_ops.py
import scanpy as sc
import anndata as ad
from ..types import Result

def pca_via_scanpy(X, n_comps: int = 50):
    tmp = ad.AnnData(X=X.copy())
    sc.pp.pca(tmp, n_comps=n_comps, use_highly_variable=False)
    return Result(outputs={"obsm.X_pca": tmp.obsm["X_pca"]},
                  state={"n_comps": n_comps})

def normalize_total_via_scanpy(X, target_sum: float = 1e4):
    tmp = ad.AnnData(X=X.copy())
    sc.pp.normalize_total(tmp, target_sum=target_sum)
    return Result(outputs={"layer.norm": tmp.X}, state={"target_sum": target_sum})

# lib/ops/scanpy_ops.py
from typing import Iterable, Optional, Sequence, Tuple, Dict, Any
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from ..types import Result

def compute_qc_metrics_scanpy(
    adata,                       # AnnData (won't be mutated)
    *,
    layer: Optional[str] = None, # which layer to use as counts (None -> .X)
    qc_vars: Optional[Sequence[str]] = None,  # names in .var that mark gene sets, e.g. ["mt"]
    percent_top: Optional[Iterable[int]] = (50, 100, 200, 500),
    by_groups: Optional[str] = None,          # obs column for per-group QC (rarely needed)
    log1p: bool = False,
    inplace: bool = True,        # we always run inplace=True on a scratch AnnData
) -> Result:
    """
    Run scanpy.pp.calculate_qc_metrics on a scratch AnnData and return a Result
    containing ONLY the newly created QC columns as DataFrames (obs & var).
    Does NOT mutate the input AnnData.
    """
    # Build scratch AnnData with the matrix we want
    X = adata.layers[layer] if layer is not None else adata.X
    tmp = ad.AnnData(X=X.copy())

    # carry over names (important for indexing/consistency)
    tmp.obs_names = adata.obs_names
    tmp.var_names = adata.var_names

    # carry over any qc_vars columns from .var if requested
    if qc_vars:
        for name in qc_vars:
            if name in adata.var:
                tmp.var[name] = adata.var[name].values
            else:
                raise KeyError(f"qc_vars entry '{name}' not found in adata.var")

    # Capture columns present before QC
    obs_before = set(tmp.obs.columns)
    var_before = set(tmp.var.columns)

    # Run scanpy (on the scratch object)
    sc.pp.calculate_qc_metrics(
        tmp,
        inplace=True,
        qc_vars=qc_vars,
        percent_top=tuple(percent_top) if percent_top is not None else None,
        log1p=log1p,
        layer=None,               # we already selected X/layer above
        # by_groups is supported in newer scanpy; if unavailable, you can remove it
        # keep a safe call:
        **({"by_groups": by_groups} if by_groups is not None else {}),
    )

    # Determine which columns were created
    obs_added = [c for c in tmp.obs.columns if c not in obs_before]
    var_added = [c for c in tmp.var.columns if c not in var_before]

    obs_qc = tmp.obs[obs_added].copy() if obs_added else pd.DataFrame(index=tmp.obs_names)
    var_qc = tmp.var[var_added].copy() if var_added else pd.DataFrame(index=tmp.var_names)

    # Package as a Result with plain pandas (no AnnData writes)
    return Result(
        outputs={
            "obs_df.qc_metrics": obs_qc,
            "var_df.qc_metrics": var_qc,
        },
        state={
            "layer": layer,
            "qc_vars": list(qc_vars) if qc_vars else None,
            "percent_top": list(percent_top) if percent_top is not None else None,
            "log1p": bool(log1p),
            "by_groups": by_groups,
        },
        metrics={
            # a couple handy aggregates (you can add more)
            "median_total_counts": float(obs_qc["total_counts"].median()) if "total_counts" in obs_qc else np.nan,
            "median_n_genes_by_counts": float(obs_qc["n_genes_by_counts"].median()) if "n_genes_by_counts" in obs_qc else np.nan,
        },
    )



from typing import Optional, Sequence, Dict, Any
import numpy as np
import anndata as ad
import scanpy as sc
from ..types import Result

def pca_from_counts_scanpy(
    adata,
    *,
    layer: Optional[str] = None,          # counts source (None -> .X)
    target_sum: float = 1e4,              # normalize_total
    log1p: bool = True,                   # log1p after normalize
    n_hvgs: int = 2000,                   # HVG count
    hvg_flavor: str = "seurat",        # "seurat", "seurat_v3", "cell_ranger"
    scale_max: Optional[float] = 10.0,    # sc.pp.scale(max_value=...)
    n_comps: int = 50,
    zero_center: bool = True,
    svd_solver: str = "arpack",
    random_state: int = 0,
) -> Result:
    """
    Functional wrapper: derives PCA embedding from counts without mutating the input.
    Returns only compact artifacts (no normalized matrix).
    """
    X = adata.layers[layer] if layer is not None else adata.X
    tmp = ad.AnnData(X=X.copy())
    tmp.obs_names = adata.obs_names
    tmp.var_names = adata.var_names

    # normalize (+ optional log1p)
    sc.pp.normalize_total(tmp, target_sum=target_sum)
    if log1p:
        sc.pp.log1p(tmp)

    # HVGs
    sc.pp.highly_variable_genes(tmp, n_top_genes=n_hvgs, flavor=hvg_flavor)
    hvg_mask = tmp.var["highly_variable"].to_numpy()
    hvg_names = tmp.var_names[hvg_mask].to_numpy()

    # subset to HVGs and (optionally) scale
    tmp = tmp[:, hvg_mask].copy()
    if scale_max is not None:
        sc.pp.scale(tmp, max_value=scale_max)

    # PCA
    sc.tl.pca(
        tmp,
        n_comps=n_comps,
        zero_center=zero_center,
        svd_solver=svd_solver,
        random_state=random_state,
        use_highly_variable=False,   # we already subset
    )

    # Collect compact outputs
    X_pca = tmp.obsm["X_pca"]                         # (cells x n_comps)
    PCs = tmp.varm["PCs"] if "PCs" in tmp.varm else None  # (HVGs x n_comps)
    pca_uns: Dict[str, Any] = dict(tmp.uns.get("pca", {}))  # explained_variance, _ratio, etc.

    outputs = {
        "obsm.X_pca": X_pca,
        "var.hvg_mask": hvg_mask,        # length = all genes
        "var.hvg_names": hvg_names,      # length = n_hvgs
        "uns.pca": pca_uns,              # small dict
    }
    if PCs is not None:
        outputs["varm.PCs"] = PCs        # loadings over HVGs only

    return Result(
        outputs=outputs,
        state={
            "layer": layer,
            "target_sum": target_sum,
            "log1p": log1p,
            "n_hvgs": n_hvgs,
            "hvg_flavor": hvg_flavor,
            "scale_max": scale_max,
            "n_comps": n_comps,
            "zero_center": zero_center,
            "svd_solver": svd_solver,
            "random_state": random_state,
        },
        metrics={
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
        },
    )


from typing import Optional, Dict, Any
import anndata as ad
import numpy as np
import scanpy as sc
from ..types import Result

def harmony_from_counts_scanpy(
    adata,
    *,
    batch_key: str,
    layer: Optional[str] = None,
    target_sum: float = 1e4,
    log1p: bool = True,
    n_hvgs: int = 2000,
    hvg_flavor: str = "seurat",
    scale_max: float | None = 10.0,
    n_comps: int = 50,
    random_state: int = 0,
) -> Result:
    """
    Normalize -> log1p -> HVG -> scale -> PCA -> Harmony.
    Returns only embedding + small metadata.
    """
    X = adata.layers[layer] if layer is not None else adata.X
    tmp = ad.AnnData(X=X.copy(), obs=adata.obs[[batch_key]].copy())
    tmp.obs_names = adata.obs_names
    tmp.var_names = adata.var_names

    sc.pp.normalize_total(tmp, target_sum=target_sum)
    if log1p:
        sc.pp.log1p(tmp)
    sc.pp.highly_variable_genes(tmp, n_top_genes=n_hvgs, flavor=hvg_flavor)
    hvg_mask = tmp.var["highly_variable"].to_numpy()
    tmp = tmp[:, hvg_mask].copy()
    if scale_max is not None:
        sc.pp.scale(tmp, max_value=scale_max)

    sc.tl.pca(tmp, n_comps=n_comps, svd_solver="arpack", random_state=random_state, use_highly_variable=False)

    # Harmony (scanpy external)
    sc.external.pp.harmony_integrate(tmp, key=batch_key, basis="X_pca")  # writes obsm["X_pca_harmony"]
    X_harmony = tmp.obsm["X_pca_harmony"]

    return Result(
        outputs={
            "obsm.X_harmony": X_harmony,
            "var.hvg_mask": hvg_mask,
        },
        state={
            "batch_key": batch_key, "layer": layer, "target_sum": target_sum,
            "log1p": log1p, "n_hvgs": n_hvgs, "hvg_flavor": hvg_flavor,
            "scale_max": scale_max, "n_comps": n_comps, "random_state": random_state,
        },
    )

# def neighbors_from_embedding_scanpy(
#     X_embed: np.ndarray,
#     *,
#     n_neighbors: int = 15,
#     metric: str = "euclidean",
#     use_rep_key: str = "X_embed",
#     random_state: int = 0,
# ) -> Result:
#     """
#     Build neighbors on a scratch AnnData with X = embedding.
#     Returns connectivities + distances (CSR matrices).
#     """
#     tmp = ad.AnnData(X=X_embed)
#     sc.pp.neighbors(tmp, n_neighbors=n_neighbors, metric=metric, use_rep=None, random_state=random_state)
#     # Scanpy stores graphs in .obsp and in .uns["neighbors"]
#     C = tmp.obsp["connectivities"]  # CSR
#     D = tmp.obsp["distances"]       # CSR
#     return Result(
#         outputs={"obsp.connectivities": C, "obsp.distances": D},
#         state={"n_neighbors": n_neighbors, "metric": metric, "use_rep": use_rep_key, "random_state": random_state},
#     )

def leiden_from_neighbors_scanpy(
    connectivities,
    *,
    resolution: float = 1.0,
    key_added: str = "leiden",
    random_state: int = 0,
) -> Result:
    """
    Leiden on a scratch AnnData by injecting the connectivities graph.
    """
    tmp = ad.AnnData(X=np.zeros((connectivities.shape[0], 1)))
    tmp.obsp["connectivities"] = connectivities
    tmp.uns["neighbors"] = {"connectivities_key": "connectivities"}
    sc.tl.leiden(tmp, resolution=resolution, key_added=key_added, random_state=random_state)
    labels = tmp.obs[key_added].astype("category").to_numpy()
    return Result(
        outputs={"obs.labels": labels},
        state={"resolution": resolution, "key_added": key_added, "random_state": random_state},
    )

# workflow/lib/ops/scanpy_ops.py

def neighbors_from_embedding_scanpy(
    X_embed: np.ndarray,
    *,
    n_neighbors: int = 15,
    metric: str = "euclidean",
    use_rep_key: str = "X_embed",
    random_state: int = 0,
) -> Result:
    """
    Build neighbors on a scratch AnnData with X = embedding.
    Returns connectivities + distances + neighbors-uns metadata.
    """
    import anndata as ad
    import scanpy as sc
    from ..types import Result

    tmp = ad.AnnData(X=X_embed)
    sc.pp.neighbors(tmp, n_neighbors=n_neighbors, metric=metric, use_rep=None, random_state=random_state)

    # graphs
    C = tmp.obsp["connectivities"]
    D = tmp.obsp["distances"]

    # # full neighbors uns (as produced by scanpy)
    # neigh_uns = dict(tmp.uns["neighbors"])  # shallow copy
    # # add a bit of context so it's clear this came from an external embedding
    # neigh_uns.setdefault("params", {}).update({"use_rep": None, "n_pcs": None})
    # neigh_uns["external_embedding_key"] = use_rep_key
    # neigh_uns["connectivities_key"] = "connectivities"
    # neigh_uns["distances_key"] = "distances"

    # full neighbors metadata exactly as produced by scanpy
    import copy
    _neigh = tmp.uns.get("neighbors", {})

    # make a deep copy so later JSON save / edits don’t mutate tmp
    neigh_uns = copy.deepcopy(_neigh)

    # (optional) ensure these keys are present & match what we saved in .obsp
    # scanpy already uses these defaults, but this keeps it explicit
    neigh_uns["connectivities_key"] = "connectivities"
    if "distances" in tmp.obsp:
        neigh_uns["distances_key"] = "distances"
    else:
        neigh_uns.pop("distances_key", None)

    return Result(
        outputs={
            "obsp.connectivities": C,
            "obsp.distances": D,
            "uns.neighbors": neigh_uns,   # <<—— new
        },
        state={"n_neighbors": n_neighbors, "metric": metric, "random_state": random_state},
    )

# def umap_from_neighbors_scanpy(
#     connectivities, distances=None, *, min_dist: float = 0.5, spread: float = 1.0, n_components: int = 2, random_state: int = 0,
# ) -> Result:
#     """
#     UMAP using precomputed neighbors.
#     """
#     n = connectivities.shape[0]
#     tmp = ad.AnnData(X=np.zeros((n, 1)))
#     tmp.obsp["connectivities"] = connectivities
#     if distances is not None:
#         tmp.obsp["distances"] = distances
#     tmp.uns["neighbors"] = {
#         "connectivities_key": "connectivities",
#         **({"distances_key": "distances"} if distances is not None else {}),
#     }
#     sc.tl.umap(tmp, min_dist=min_dist, spread=spread, n_components=n_components, random_state=random_state)
#     X_umap = tmp.obsm["X_umap"]
#     return Result(
#         outputs={"obsm.X_umap": X_umap},
#         state={"min_dist": min_dist, "spread": spread, "n_components": n_components, "random_state": random_state},
#     )

# def umap_from_neighbors_scanpy(
#     connectivities, distances=None, *,
#     min_dist: float = 0.5, spread: float = 1.0, n_components: int = 2,
#     random_state: int = 0,
# ) -> Result:
#     """
#     UMAP using precomputed neighbors graphs (connectivities[, distances]).
#     We must also provide a minimal neighbors params dict expected by Scanpy.
#     """
#     import anndata as ad
#     import numpy as np
#     import scanpy as sc
#     from ..types import Result

#     n = connectivities.shape[0]
#     tmp = ad.AnnData(X=np.zeros((n, 1)))

#     tmp.obsp["connectivities"] = connectivities
#     if distances is not None:
#         tmp.obsp["distances"] = distances

#     neigh_uns = {
#         "connectivities_key": "connectivities",
#         **({"distances_key": "distances"} if distances is not None else {}),
#         # Minimal params Scanpy's UMAP reads:
#         "params": {
#             "use_rep": None,     # we’re not using a data rep, graphs are provided
#             "n_pcs": None,
#             # (optional extras; harmless but informative)
#             "method": "umap",
#             "n_neighbors": None,
#             "metric": "precomputed",
#         },
#     }
#     tmp.uns["neighbors"] = neigh_uns

#     sc.tl.umap(tmp, min_dist=min_dist, spread=spread,
#                n_components=n_components, random_state=random_state)
#     X_umap = tmp.obsm["X_umap"]

#     return Result(
#         outputs={"obsm.X_umap": X_umap},
#         state={"min_dist": min_dist, "spread": spread,
#                "n_components": n_components, "random_state": random_state},
#     )

def umap_from_neighbors_scanpy(
    connectivities, distances=None, *,
    min_dist: float = 0.5, spread: float = 1.0, n_components: int = 2,
    random_state: int = 0, neighbors_uns: dict | None = None
) -> Result:
    import anndata as ad, numpy as np, scanpy as sc
    from ..types import Result

    n = connectivities.shape[0]
    tmp = ad.AnnData(X=np.zeros((n, 1)))
    tmp.obsp["connectivities"] = connectivities
    if distances is not None:
        tmp.obsp["distances"] = distances

    # if caller didn’t provide metadata, create a minimal, valid one
    if neighbors_uns is None:
        neighbors_uns = {
            "connectivities_key": "connectivities",
            **({"distances_key": "distances"} if distances is not None else {}),
            "params": {
                "method": "umap",
                "use_rep": None,
                "n_pcs": None,
                "metric": "precomputed",
            },
        }
    else:
        # ensure keys match obsp keys and required params exist
        neighbors_uns = dict(neighbors_uns)  # shallow copy is fine
        neighbors_uns["connectivities_key"] = "connectivities"
        if distances is not None:
            neighbors_uns["distances_key"] = "distances"
        params = dict(neighbors_uns.get("params", {}))
        params.setdefault("method", "umap")
        params.setdefault("use_rep", None)
        params.setdefault("n_pcs", None)
        neighbors_uns["params"] = params

    tmp.uns["neighbors"] = neighbors_uns

    sc.tl.umap(tmp, min_dist=min_dist, spread=spread,
               n_components=n_components, random_state=random_state)
    return Result(outputs={"obsm.X_umap": tmp.obsm["X_umap"]},
                  state={"min_dist": min_dist, "spread": spread,
                         "n_components": n_components, "random_state": random_state})