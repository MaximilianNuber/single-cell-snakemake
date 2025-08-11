from __future__ import annotations
import anndata as ad
import scanpy as sc
from ...types import Result

class HarmonyScanpy:
    """
    counts -> normalize/log/HVG/scale -> PCA -> Harmony (scanpy.external)
    Returns obsm.X_harmony
    """
    def __init__(
        self,
        batch_key: str,
        n_comps: int = 50,
        target_sum: float = 1e4,
        log1p: bool = True,
        n_hvgs: int = 2000,
        hvg_flavor: str = "seurat_v3",
        scale_max: float | None = 10.0,
        random_state: int = 0,
    ):
        self.cfg = dict(
            batch_key=batch_key, n_comps=n_comps, target_sum=target_sum,
            log1p=log1p, n_hvgs=n_hvgs, hvg_flavor=hvg_flavor,
            scale_max=scale_max, random_state=random_state,
        )

    def compute(self, X_counts, obs, var) -> Result:
        try:
            import harmonypy  # noqa: F401
        except ImportError as e:
            raise ImportError("Harmony requires `harmonypy` (pip install harmonypy).") from e

        tmp = ad.AnnData(X=X_counts, obs=obs[[self.cfg["batch_key"]]].copy(), var=var.copy())

        sc.pp.normalize_total(tmp, target_sum=self.cfg["target_sum"])
        if self.cfg["log1p"]:
            sc.pp.log1p(tmp)
        sc.pp.highly_variable_genes(tmp, n_top_genes=self.cfg["n_hvgs"], flavor=self.cfg["hvg_flavor"])
        tmp = tmp[:, tmp.var["highly_variable"].to_numpy()].copy()
        if self.cfg["scale_max"] is not None:
            sc.pp.scale(tmp, max_value=self.cfg["scale_max"])

        sc.tl.pca(tmp, n_comps=self.cfg["n_comps"], random_state=self.cfg["random_state"], use_highly_variable=False)
        sc.external.pp.harmony_integrate(tmp, key=self.cfg["batch_key"], basis="X_pca")  # => obsm["X_pca_harmony"]

        return Result(outputs={"obsm.X_harmony": tmp.obsm["X_pca_harmony"]}, state=self.cfg)