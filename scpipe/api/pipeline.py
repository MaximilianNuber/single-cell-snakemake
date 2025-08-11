from __future__ import annotations
import numpy as np
from typing import Any
from ..protocols import EmbeddingComputer, NeighborGraphBuilder, ManifoldEmbedder

def run_pipeline(
    X_counts: Any, obs: Any, var: Any,
    emb: EmbeddingComputer,
    nn: NeighborGraphBuilder,
    manifold: ManifoldEmbedder,
):
    assert isinstance(emb, EmbeddingComputer), "emb must implement .compute(...) -> Result"
    Z = emb.compute(X_counts, obs, var).outputs["obsm.X_pca"]  # or X_harmony/scvi

    assert isinstance(nn, NeighborGraphBuilder), "nn must implement .build(...) -> Result"
    r_nn = nn.build(Z)
    C, D = r_nn.outputs["obsp.connectivities"], r_nn.outputs["obsp.distances"]

    assert isinstance(manifold, ManifoldEmbedder), "manifold must implement .embed(...) -> Result"
    U = manifold.embed(C, distances=D).outputs["obsm.X_umap"]

    return {"embedding": Z, "connectivities": C, "distances": D, "umap": U}