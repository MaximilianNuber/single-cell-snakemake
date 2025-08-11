# workflow/scripts/run_scdblfinder.py
import pandas as pd
import anndata as ad
from scpipe.interop.r_qc import run_scdblfinder

h5ad_in = snakemake.input["h5ad"]
tsv_out = snakemake.output["tsv"]
batch_key = snakemake.params["batch_key"]

adata = ad.read_h5ad(h5ad_in)
adata = run_scdblfinder(adata, batch_key=batch_key)

df = pd.DataFrame({
    "barcode": adata.obs_names.astype(str),
    "is_doublet": adata.obs["is_doublet"].astype(bool).to_numpy(),
    "dbl_score": adata.obs.get("dbl_score", pd.Series(index=adata.obs_names, dtype=float)).to_numpy(),
})
df.to_csv(tsv_out, sep="\t", index=False)
