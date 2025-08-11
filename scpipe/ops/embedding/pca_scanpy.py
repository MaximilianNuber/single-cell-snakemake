from __future__ import annotations
import scanpy as sc
import anndata as ad
import numpy as np
from ...types import Result  # adjust if your relative path differs

class PCAScanpy:
    """
    counts -> (normalize -> log1p) -> HVG -> (scale) -> PCA
    Returns only the embedding (obsm.X_pca) and small metadata.
    """
    def __init__(
        self,
        n_comps: int = 50,
        target_sum: float = 1e4,
        log1p: bool = True,
        n_hvgs: int = 2000,
        hvg_flavor: str = "seurat_v3",
        scale_max: float | None = 10.0,
        random_state: int = 0,
    ):
        self.cfg = dict(
            n_comps=n_comps, target_sum=target_sum, log1p=log1p,
            n_hvgs=n_hvgs, hvg_flavor=hvg_flavor, scale_max=scale_max,
            random_state=random_state,
        )

    def compute(self, X_counts, obs, var) -> Result:
        tmp = ad.AnnData(X=X_counts, obs=obs.copy(), var=var.copy())

        sc.pp.normalize_total(tmp, target_sum=self.cfg["target_sum"])
        if self.cfg["log1p"]:
            sc.pp.log1p(tmp)

        sc.pp.highly_variable_genes(tmp, n_top_genes=self.cfg["n_hvgs"], flavor=self.cfg["hvg_flavor"])
        tmp = tmp[:, tmp.var["highly_variable"].to_numpy()].copy()

        if self.cfg["scale_max"] is not None:
            sc.pp.scale(tmp, max_value=self.cfg["scale_max"])

        sc.tl.pca(tmp, n_comps=self.cfg["n_comps"], random_state=self.cfg["random_state"], use_highly_variable=False)

        return Result(
            outputs={
                "obsm.X_pca": tmp.obsm["X_pca"],
                "uns.pca": dict(tmp.uns.get("pca", {})),
                "var.hvg_names": tmp.var_names.to_numpy(),
            },
            state=self.cfg,
        )