from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad
from scpipe.ops.artifacts import save_embedding, save_json
from scpipe.ops.integration_ops import (
    latent_from_pca_identity,
    latent_from_pca_harmony,
    latent_from_counts_scvi,
)

run_cfg = dict(snakemake.params["run_cfg"])
method = (run_cfg.get("method") or "none").lower()
batch_key = run_cfg.get("batch_key", "sample")

A = ad.read_h5ad(snakemake.input["concat"])
X_pca = None
if "pca" in snakemake.input and snakemake.input["pca"]:
    X_pca = np.load(snakemake.input["pca"])

# choose path
if method in {"none", "unintegrated"}:
    res = latent_from_pca_identity(X_pca)
elif method == "bbknn":
    # embedding is PCA, neighbors backend switches later
    res = latent_from_pca_identity(X_pca)
    # annotate in state for downstream
    res.state = {**res.state, "neighbors_backend": "bbknn", "batch_key": batch_key}
elif method == "harmony":
    res = latent_from_pca_harmony(X_pca, obs=A.obs, batch_key=batch_key,
                                  random_state=int(run_cfg.get("random_state", 0)))
elif method == "scvi":
    res = latent_from_counts_scvi(A, batch_key=batch_key,
                                  n_latent=int(run_cfg.get("n_latent", 30)),
                                  max_epochs=int(run_cfg.get("max_epochs", 200)),
                                  seed=int(run_cfg.get("seed", 0)))
else:
    raise ValueError(f"Unknown integration method: {method}")

# write artifacts
out_dir = Path(snakemake.output["manifest"]).parent
out_dir.mkdir(parents=True, exist_ok=True)
paths = save_embedding(res.outputs["obsm.X_latent"], A.obs_names, out_dir, key="X_latent")

manifest = {
    "method": method,
    "embedding_key": res.outputs.get("embedding_key", "X_latent"),
    "neighbors_backend": res.state.get("neighbors_backend", "standard"),
    "batch_key": batch_key,
    "params": {k: v for k, v in res.state.items()
               if k not in {"method", "batch_key", "neighbors_backend"}},
}
save_json(manifest, snakemake.output["manifest"])
