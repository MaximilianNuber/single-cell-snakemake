# lib/api/graph.py
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from scpipe.ops.scanpy_ops import neighbors_from_embedding_scanpy, leiden_from_neighbors_scanpy, umap_from_neighbors_scanpy
from ..ops.artifacts import save_sparse_csr, save_labels, save_embedding, save_json

def _load_embedding(path_npy: str, path_index: str, target_obs_names) -> np.ndarray:
    X = np.load(path_npy)
    idx = pd.read_csv(path_index)["obs_names"].astype(str).to_numpy()
    pos = {n: i for i, n in enumerate(idx)}
    order = np.array([pos[n] for n in target_obs_names], dtype=int)
    return X[order]

# def compute_and_save_neighbors_from_embedding(
#     obs_names,
#     embedding_paths: dict,          # {"path": "...npy", "index": "...csv"}
#     out_dir,
#     *,
#     n_neighbors: int = 15,
#     metric: str = "euclidean",
# ) -> dict:
#     X = _load_embedding(embedding_paths["path"], embedding_paths["index"], obs_names)
#     res = neighbors_from_embedding_scanpy(X, n_neighbors=n_neighbors, metric=metric, use_rep_key="external")
#     out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
#     paths = {
#         "obsp.connectivities": save_sparse_csr(res.outputs["obsp.connectivities"], out / "connectivities.npz"),
#         "obsp.distances":      save_sparse_csr(res.outputs["obsp.distances"], out / "distances.npz"),
#     }
#     save_json({"params": res.state}, out / "neighbors.json")
#     return {"graphs": paths}

def compute_and_save_neighbors_from_embedding(
    obs_names,
    embedding_paths: dict,          # {"path": "...npy", "index": "...csv"}
    out_dir,
    *,
    n_neighbors: int = 15,
    metric: str = "euclidean",
) -> dict:
    X = _load_embedding(embedding_paths["path"], embedding_paths["index"], obs_names)
    res = neighbors_from_embedding_scanpy(X, n_neighbors=n_neighbors, metric=metric, use_rep_key="external")

    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    paths = {
        "obsp.connectivities": save_sparse_csr(res.outputs["obsp.connectivities"], out / "connectivities.npz"),
        "obsp.distances":      save_sparse_csr(res.outputs["obsp.distances"], out / "distances.npz"),
    }

    # save the neighbors metadata exactly as Scanpy produced it
    neighbors_json = save_json(res.outputs["uns.neighbors"], out / "neighbors.json")

    return {"graphs": paths, "uns.neighbors": neighbors_json}

def compute_and_save_leiden_from_neighbors(
    obs_names,
    graph_paths: dict,              # {"connectivities": "...npz"}
    out_csv,
    *,
    resolution: float = 1.0,
    key: str = "leiden",
) -> dict:
    from scipy.sparse import load_npz
    C = load_npz(graph_paths["connectivities"])
    res = leiden_from_neighbors_scanpy(C, resolution=resolution, key_added=key)
    path = save_labels(res.outputs["obs.labels"], obs_names, out_csv, key=key)
    return {"obs.labels": path}

# def compute_and_save_umap_from_neighbors(
#     obs_names,
#     graph_paths: dict,              # {"connectivities": "...npz", "distances": "...npz" (optional)}
#     out_dir,
#     *,
#     min_dist: float = 0.5,
#     spread: float = 1.0,
#     n_components: int = 2,
# ) -> dict:
#     from scipy.sparse import load_npz
#     C = load_npz(graph_paths["connectivities"])
#     D = load_npz(graph_paths["distances"]) if "distances" in graph_paths else None
#     res = umap_from_neighbors_scanpy(C, distances=D, min_dist=min_dist, spread=spread, n_components=n_components)
#     paths = {}
#     paths["obsm.X_umap"] = save_embedding(res.outputs["obsm.X_umap"], obs_names, out_dir, key="X_umap")
#     return paths

def compute_and_save_umap_from_neighbors(
    obs_names,
    graph_paths: dict,              # {"connectivities": "...npz", "distances": "...npz"}
    out_dir,
    *,
    neighbors_json: str | None = None,   # <<—— pass the JSON path we saved
    min_dist: float = 0.5,
    spread: float = 1.0,
    n_components: int = 2,
) -> dict:
    from scipy.sparse import load_npz
    import json

    C = load_npz(graph_paths["connectivities"])
    D = load_npz(graph_paths["distances"]) if "distances" in graph_paths else None

    # Load neighbors uns metadata if provided; otherwise fallback to a minimal dict
    if neighbors_json:
        with open(neighbors_json) as f:
            neigh_uns = json.load(f)
    else:
        neigh_uns = {
            "connectivities_key": "connectivities",
            **({"distances_key": "distances"} if D is not None else {}),
            "params": {"use_rep": None, "n_pcs": None, "metric": "precomputed"},
        }

    # Hand the full neighbors-uns to ops.umap
    res = umap_from_neighbors_scanpy(
        C, distances=D, min_dist=min_dist, spread=spread, n_components=n_components,
        # we'll just set it inside using the same keys:
    )

    # Save embedding
    paths = {}
    paths["obsm.X_umap"] = save_embedding(res.outputs["obsm.X_umap"], obs_names, out_dir, key="X_umap")
    return paths