# scpipe/api/embedding.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, Sequence
from ..types import Result
from ..ops.artifacts import save_embedding, save_loadings, save_json

def _choose_emb_key(outputs: Dict[str, Any]) -> str:
    for k in ("obsm.X_pca", "obsm.X_scvi", "obsm.X_pca_harmony"):
        if k in outputs:
            return k
    raise ValueError("Embedding Result lacks a known embedding key (obsm.X_pca / X_scvi / X_pca_harmony).")

def save_embedding_result(res: Result, obs_names: Sequence[str], out_dir: str | Path) -> Dict[str, Any]:
    if res.kind != "embedding":
        raise ValueError(f"save_embedding_result expects Result(kind='embedding'), got {res.kind}")
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)

    emb_key = _choose_emb_key(res.outputs)
    short = emb_key.split(".")[-1]  # e.g., X_pca
    paths: Dict[str, Any] = {}
    paths[emb_key] = save_embedding(res.outputs[emb_key], obs_names, out, key=short)

    if res.outputs.get("varm.PCs") is not None and res.outputs.get("var.hvg_names") is not None:
        paths["varm.PCs"] = save_loadings(res.outputs["varm.PCs"], res.outputs["var.hvg_names"], out, key="PCs")

    manifest = {
        "kind": res.kind,
        "params": res.state,
        "metrics": res.metrics,
        "uns": {"pca": res.outputs.get("uns.pca")},
    }
    paths["manifest"] = {"json": save_json(manifest, out / "manifest.json")}
    return paths

# from __future__ import annotations
# # lib/api/embedding.py

# from typing import Optional
# from pathlib import Path
# from ..ops.scanpy_ops import pca_from_counts_scanpy
# from ..ops.artifacts import save_embedding, save_loadings, save_json

# def compute_and_save_pca(
#     adata,
#     out_dir: str | Path,
#     *,
#     layer: Optional[str] = None,
#     target_sum: float = 1e4,
#     log1p: bool = True,
#     n_hvgs: int = 2000,
#     hvg_flavor: str = "seurat_v3",
#     scale_max: float | None = 10.0,
#     n_comps: int = 50,
#     zero_center: bool = True,
#     svd_solver: str = "arpack",
#     random_state: int = 0,
#     save_loadings_too: bool = False,
# ) -> dict:
#     """
#     Derives PCA purely from counts and saves:
#       - obsm.X_pca.npy + index CSV (obs_names)
#       - pca.json with small metadata (explained variance, params, HVG names)
#       - (optional) PCs.npy + genes CSV over HVGs
#     """
#     res = pca_from_counts_scanpy(
#         adata,
#         layer=layer,
#         target_sum=target_sum,
#         log1p=log1p,
#         n_hvgs=n_hvgs,
#         hvg_flavor=hvg_flavor,
#         scale_max=scale_max,
#         n_comps=n_comps,
#         zero_center=zero_center,
#         svd_solver=svd_solver,
#         random_state=random_state,
#     )

#     out_dir = Path(out_dir)
#     paths = {}
#     # Save embedding
#     paths["obsm.X_pca"] = save_embedding(
#         res.outputs["obsm.X_pca"], adata.obs_names, out_dir, key="X_pca"
#     )

#     # Save small PCA metadata JSON
#     pca_meta = {
#         "explained_variance": res.outputs.get("uns.pca", {}).get("variance", None),
#         "explained_variance_ratio": res.outputs.get("uns.pca", {}).get("variance_ratio", None),
#         "hvg_names": res.outputs.get("var.hvg_names", []).tolist() if "var.hvg_names" in res.outputs else [],
#         "params": res.state,
#     }
#     paths["uns.pca"] = {"json": save_json(pca_meta, out_dir / "pca.json")}

#     # Optionally save loadings on HVGs
#     if save_loadings_too and "varm.PCs" in res.outputs:
#         paths["varm.PCs"] = save_loadings(
#             res.outputs["varm.PCs"], res.outputs["var.hvg_names"], out_dir, key="PCs"
#         )

#     return paths


# # lib/api/embedding.py
# from pathlib import Path
# from typing import Optional
# from ..ops.scanpy_ops import harmony_from_counts_scanpy
# from ..ops.scvi_ops import scvi_from_counts
# from ..ops.artifacts import save_embedding, save_json

# def compute_and_save_harmony(adata, out_dir, *, batch_key, **kwargs) -> dict:
#     res = harmony_from_counts_scanpy(adata, batch_key=batch_key, **kwargs)
#     out = Path(out_dir)
#     paths = {}
#     paths["obsm.X_harmony"] = save_embedding(res.outputs["obsm.X_harmony"], adata.obs_names, out, key="X_harmony")
#     paths["uns.harmony"] = {"json": save_json({"params": res.state}, out / "harmony.json")}
#     return paths

# def compute_and_save_scvi(adata, out_dir, *, batch_key, **kwargs) -> dict:
#     res = scvi_from_counts(adata, batch_key=batch_key, **kwargs)
#     out = Path(out_dir)
#     paths = {}
#     paths["obsm.X_scvi"] = save_embedding(res.outputs["obsm.X_scvi"], adata.obs_names, out, key="X_scvi")
#     paths["uns.scvi"] = {"json": save_json({"params": res.state}, out / "scvi.json")}
#     return paths
