# scpipe/ops/integration_ops.py
from __future__ import annotations
from typing import Any, Dict, Mapping
import numpy as np
import pandas as pd
import anndata as ad
from ..types import Result

def latent_from_pca_identity(X_pca: np.ndarray, *, var_names=None, prep_state: Dict[str, Any] | None = None) -> Result:
    """No integration: pass PCA through as latent."""
    return Result(
        kind="embedding",
        outputs={"obsm.X_latent": X_pca, "embedding_key": "X_pca"},
        state={"schema": "latent/identity@v1", "neighbors_backend": "standard", **(prep_state or {})},
        metrics={},
    )

def latent_from_pca_harmony(
    X_pca: np.ndarray,
    *,
    obs: Mapping[str, Any],
    batch_key: str = "sample",
    random_state: int = 0,
    prep_state: Dict[str, Any] | None = None,
) -> Result:
    """Harmony on PCA; returns X_latent = X_pca_harmony."""
    import scanpy as sc
    import scanpy.external as sce
    A = ad.AnnData(X=np.zeros((len(obs), 1)))  # throwaway; we only store embeddings in obsm
    A.obs = pd.DataFrame(obs).copy()
    A.obsm["X_pca"] = np.asarray(X_pca)
    sce.pp.harmony_integrate(A, key=batch_key, random_state=random_state)
    Xh = A.obsm["X_pca_harmony"]
    return Result(
        kind="embedding",
        outputs={"obsm.X_latent": Xh, "embedding_key": "X_pca_harmony"},
        state={"schema": "latent/harmony@v1", "batch_key": batch_key, "neighbors_backend": "standard",
               "random_state": random_state, **(prep_state or {})},
        metrics={},
    )

def latent_from_counts_scvi(
    adata: ad.AnnData,
    *,
    batch_key: str = "sample",
    n_latent: int = 30,
    max_epochs: int = 200,
    seed: int = 0,
    extra_state: Dict[str, Any] | None = None,
) -> Result:
    """scVI on counts; returns X_scvi as X_latent."""
    from .scvi_ops import scvi_from_counts  # your existing op returning Result with obsm.X_scvi
    r = scvi_from_counts(adata, batch_key=batch_key, n_latent=n_latent, max_epochs=max_epochs, seed=seed)
    Z = r.outputs["obsm.X_scvi"]
    return Result(
        kind="embedding",
        outputs={"obsm.X_latent": Z, "embedding_key": "X_scvi"},
        state={"schema": "latent/scvi@v1", "batch_key": batch_key, "n_latent": n_latent,
               "max_epochs": max_epochs, "seed": seed, "neighbors_backend": "standard", **(extra_state or {})},
        metrics=r.metrics,
    )
