from __future__ import annotations
import anndata as ad
from ...types import Result

class SCVIEmbedding:
    """
    Fit scVI on counts and return obsm.X_scvi (latent).
    """
    def __init__(self, batch_key: str, n_latent: int = 30, max_epochs: int = 200, seed: int = 0):
        self.cfg = dict(batch_key=batch_key, n_latent=n_latent, max_epochs=max_epochs, seed=seed)

    def compute(self, X_counts, obs, var) -> Result:
        try:
            import scvi
        except ImportError as e:
            raise ImportError("scvi-tools not installed. Run in env with `pip install scvi-tools`.") from e

        scvi.settings.seed = self.cfg["seed"]
        tmp = ad.AnnData(X=X_counts, obs=obs[[self.cfg["batch_key"]]].copy(), var=var.copy())
        scvi.model.SCVI.setup_anndata(tmp, batch_key=self.cfg["batch_key"])
        model = scvi.model.SCVI(tmp, n_latent=self.cfg["n_latent"])
        model.train(max_epochs=self.cfg["max_epochs"], check_val_every_n_epoch=None, enable_progress_bar=False)
        Z = model.get_latent_representation()
        return Result(outputs={"obsm.X_scvi": Z}, state=self.cfg)