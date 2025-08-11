# lib/ops/scvi_ops.py
from __future__ import annotations
from typing import Optional, Dict, Any
import numpy as np
import anndata as ad
from ..types import Result

def scvi_from_counts(
    adata,
    *,
    batch_key: str,
    layer: Optional[str] = None,        # counts source
    n_latent: int = 30,
    max_epochs: int = 200,
    seed: int = 0,
) -> Result:
    """
    Fit SCVI on counts and return the latent embedding only.
    Does not persist the model unless you decide to.
    """

    try:
        import scvi
    except ImportError as e:
        raise ImportError(
            "scvi-tools is not installed. Install it (e.g. `pip install scvi-tools`) "
            "or run this step in an environment that has it."
        ) from e

    scvi.settings.seed = seed
    X = adata.layers[layer] if layer is not None else adata.X
    tmp = ad.AnnData(X=X.copy(), obs=adata.obs[[batch_key]].copy(), var=adata.var.copy())
    tmp.obs_names = adata.obs_names
    tmp.var_names = adata.var_names

    scvi.model.SCVI.setup_anndata(tmp, layer=None, batch_key=batch_key)  # tmp.X is counts
    model = scvi.model.SCVI(tmp, n_latent=n_latent)
    model.train(max_epochs=max_epochs, check_val_every_n_epoch=None, enable_progress_bar=False)
    Z = model.get_latent_representation()  # (cells x n_latent)

    return Result(
        outputs={"obsm.X_scvi": Z},
        state={"batch_key": batch_key, "layer": layer, "n_latent": n_latent, "max_epochs": max_epochs, "seed": seed},
    )
