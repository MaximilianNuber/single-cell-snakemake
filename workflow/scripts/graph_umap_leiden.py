from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scanpy.external as sce
from scipy.sparse import save_npz
from scpipe.ops.artifacts import save_labels, save_json

A = ad.read_h5ad(snakemake.input["concat"])

X = np.load(snakemake.input["latent"])
idx = pd.read_csv(snakemake.input["idx"])["obs_names"].astype(str).to_numpy()
pos = {n: i for i, n in enumerate(idx)}
X = X[[pos[n] for n in A.obs_names]]  # align to AnnData row order
A.obsm["X_latent"] = X

# integration manifest
with open(snakemake.input["manifest"]) as f:
    integ_meta = json.load(f)
method = integ_meta.get("method", "none")
neighbors_backend = integ_meta.get("neighbors_backend", "standard")
batch_key = integ_meta.get("batch_key", "sample")

cfg = dict(snakemake.params.get("run_cfg", {}))
k = int(cfg.get("neighbors_k", 15))
min_dist = float(cfg.get("umap_min_dist", 0.30))
resolution = float(cfg.get("cluster_resolution", 1.0))

# neighbors
if neighbors_backend == "bbknn":
    n_batches = int(A.obs[batch_key].nunique()) or 1
    within = max(3, k // n_batches)
    try:
        sce.pp.bbknn(A, batch_key=batch_key, neighbors_within_batch=within, use_rep="X_latent")
    except TypeError:
        sce.pp.bbknn(A, batch_key=batch_key, neighbors_within_batch=within)
else:
    sc.pp.neighbors(A, use_rep="X_latent", n_neighbors=k)

# UMAP + Leiden
sc.tl.umap(A, min_dist=min_dist)
sc.tl.leiden(A, resolution=resolution, key_added="leiden")

# save artifacts
Path(snakemake.output["knn"]).parent.mkdir(parents=True, exist_ok=True)
save_npz(snakemake.output["knn"], A.obsp["connectivities"])

pd.DataFrame(
    {"u1": A.obsm["X_umap"][:, 0], "u2": A.obsm["X_umap"][:, 1]},
    index=A.obs_names,
).to_csv(snakemake.output["umap"])

save_labels(A.obs["leiden"].to_numpy(), A.obs_names, snakemake.output["leiden"], key="leiden")

save_json({
    "method": method,
    "neighbors_backend": neighbors_backend,
    "batch_key": batch_key,
    "neighbors_k": k,
    "umap_min_dist": min_dist,
    "cluster_resolution": resolution
}, snakemake.output["manifest"])
