# scpipe/ops/pca_ops.py
from __future__ import annotations
from functools import singledispatch
from typing import Any, Mapping, Tuple, Dict
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from ..types import Result

# ------------------------------------------------------------
# 0) Helpers
# ------------------------------------------------------------
def _ensure_df(x, n: int, axis: str) -> pd.DataFrame:
    """Coerce mapping-like into a DataFrame of length n."""
    if isinstance(x, pd.DataFrame):
        df = x.copy()
    else:
        df = pd.DataFrame(x)
    if len(df) != n:
        raise ValueError(f"{axis} length mismatch: got {len(df)}, expected {n}.")
    return df

# ------------------------------------------------------------
# 1) Shared PREP: normalize → log1p → HVG → (optional) scale
#    Returns (prepped_matrix, var_names, prep_state)
# ------------------------------------------------------------
@singledispatch
def prep_counts_for_pca(
    X_counts: Any,
    *,
    obs: Mapping,
    var: Mapping,
    target_sum: float = 1e4,
    log1p: bool = True,
    n_hvgs: int = 2000,
    hvg_flavor: str = "seurat_v3",
    scale_max: float | None = 10.0,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """Default overload: triplet interface (X, obs, var)."""
    # Build an AnnData to reuse Scanpy's prep ops
    X = X_counts
    n_obs = X.shape[0]
    n_var = X.shape[1]
    obs_df = _ensure_df(obs, n_obs, "obs")
    var_df = _ensure_df(var, n_var, "var")
    adata = ad.AnnData(X=X, obs=obs_df, var=var_df)
    return prep_counts_for_pca(adata,
                               target_sum=target_sum,
                               log1p=log1p,
                               n_hvgs=n_hvgs,
                               hvg_flavor=hvg_flavor,
                               scale_max=scale_max)

@prep_counts_for_pca.register
def _(adata: ad.AnnData,
      *,
      target_sum: float = 1e4,
      log1p: bool = True,
      n_hvgs: int = 2000,
      hvg_flavor: str = "seurat_v3",
      scale_max: float | None = 10.0,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """AnnData overload."""
    tmp = adata.copy()
    sc.pp.normalize_total(tmp, target_sum=target_sum)
    if log1p:
        sc.pp.log1p(tmp)

    sc.pp.highly_variable_genes(tmp, n_top_genes=n_hvgs, flavor=hvg_flavor)
    tmp = tmp[:, tmp.var["highly_variable"].to_numpy()].copy()

    if scale_max is not None:
        sc.pp.scale(tmp, max_value=scale_max)  # will densify if needed

    # Return a plain matrix + var_names + prep metadata
    X_prepped = tmp.X.toarray() if hasattr(tmp.X, "toarray") else np.asarray(tmp.X)
    var_names = tmp.var_names.to_numpy()
    prep_state = {
        "schema": "prep/pca@v1",
        "target_sum": target_sum,
        "log1p": log1p,
        "n_hvgs": n_hvgs,
        "hvg_flavor": hvg_flavor,
        "scale_max": scale_max,
    }
    return X_prepped, var_names, prep_state


# ------------------------------------------------------------
# 2) PCA backends: both accept a 2D matrix + var_names
# ------------------------------------------------------------
def _pca_impl_scanpy(
    X_prepped: np.ndarray,
    var_names: np.ndarray,
    *,
    n_comps: int = 50,
    zero_center: bool = True,
    svd_solver: str = "arpack",
    random_state: int = 0,
    prep_state: Dict[str, Any] | None = None,
) -> Result:
    """Run PCA via scanpy.tl.pca on a throwaway AnnData."""
    tmp = ad.AnnData(X=X_prepped, var=pd.DataFrame(index=pd.Index(var_names, name="var_names")))
    sc.tl.pca(
        tmp,
        n_comps=n_comps,
        zero_center=zero_center,
        svd_solver=svd_solver,
        random_state=random_state,
        use_highly_variable=False,
    )
    pca_uns = dict(tmp.uns.get("pca", {}))
    outputs = {
        "obsm.X_pca": tmp.obsm["X_pca"],  # (cells x PCs)
        "uns.pca": {
            "variance": pca_uns.get("variance"),
            "variance_ratio": pca_uns.get("variance_ratio"),
        },
        "var.hvg_names": var_names,
        "varm.PCs": tmp.varm.get("PCs"),  # (genes x PCs)
    }
    state = {
        "schema": "embedding/pca_scanpy@v1",
        "n_comps": n_comps,
        "zero_center": zero_center,
        "svd_solver": svd_solver,
        "random_state": random_state,
        **(prep_state or {}),
    }
    vr = outputs["uns.pca"]["variance_ratio"]
    metrics = {"pca_pc1_var_ratio": float(vr[0])} if vr is not None and len(vr) else {}
    return Result(kind="embedding", outputs=outputs, state=state, metrics=metrics)

def _pca_impl_sklearn(
    X_prepped: np.ndarray,
    var_names: np.ndarray,
    *,
    n_comps: int = 50,
    svd_solver: str = "full",
    random_state: int = 0,
    prep_state: Dict[str, Any] | None = None,
) -> Result:
    """Run PCA via scikit-learn on a dense matrix."""
    from sklearn.decomposition import PCA
    X = np.asarray(X_prepped)
    pca = PCA(n_components=n_comps, svd_solver=svd_solver, random_state=random_state)
    X_pca = pca.fit_transform(X)        # (cells x PCs)
    PCs = pca.components_.T             # (genes x PCs), to mirror scanpy's varm["PCs"]

    outputs = {
        "obsm.X_pca": X_pca,
        "uns.pca": {
            "variance": pca.explained_variance_,
            "variance_ratio": pca.explained_variance_ratio_,
        },
        "var.hvg_names": var_names,
        "varm.PCs": PCs,
    }
    state = {
        "schema": "embedding/pca_sklearn@v1",
        "n_comps": n_comps,
        "svd_solver": svd_solver,
        "random_state": random_state,
        **(prep_state or {}),
    }
    metrics = {"pca_pc1_var_ratio": float(outputs["uns.pca"]["variance_ratio"][0])}
    return Result(kind="embedding", outputs=outputs, state=state, metrics=metrics)


# ------------------------------------------------------------
# 3) Public APIs: accept AnnData OR (X, obs, var), reuse prep → backend
# ------------------------------------------------------------
@singledispatch
def pca_from_counts_scanpy(
    X_counts: Any,
    *,
    obs: Mapping,
    var: Mapping,
    **cfg,
) -> Result:
    X, vnames, pst = prep_counts_for_pca(X_counts, obs=obs, var=var, **cfg)
    return _pca_impl_scanpy(
        X, vnames,
        n_comps=cfg.get("n_comps", 50),
        zero_center=cfg.get("zero_center", True),
        svd_solver=cfg.get("svd_solver", "arpack"),
        random_state=cfg.get("random_state", 0),
        prep_state=pst,
    )

@pca_from_counts_scanpy.register
def _(adata: ad.AnnData, **cfg) -> Result:
    X, vnames, pst = prep_counts_for_pca(adata, **cfg)
    return _pca_impl_scanpy(
        X, vnames,
        n_comps=cfg.get("n_comps", 50),
        zero_center=cfg.get("zero_center", True),
        svd_solver=cfg.get("svd_solver", "arpack"),
        random_state=cfg.get("random_state", 0),
        prep_state=pst,
    )

@singledispatch
def pca_from_counts_sklearn(
    X_counts: Any,
    *,
    obs: Mapping,
    var: Mapping,
    **cfg,
) -> Result:
    X, vnames, pst = prep_counts_for_pca(X_counts, obs=obs, var=var, **cfg)
    return _pca_impl_sklearn(
        X, vnames,
        n_comps=cfg.get("n_comps", 50),
        svd_solver=cfg.get("svd_solver", "full"),
        random_state=cfg.get("random_state", 0),
        prep_state=pst,
    )

@pca_from_counts_sklearn.register
def _(adata: ad.AnnData, **cfg) -> Result:
    X, vnames, pst = prep_counts_for_pca(adata, **cfg)
    return _pca_impl_sklearn(
        X, vnames,
        n_comps=cfg.get("n_comps", 50),
        svd_solver=cfg.get("svd_solver", "full"),
        random_state=cfg.get("random_state", 0),
        prep_state=pst,
    )
