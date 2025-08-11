# workflow/scripts/run_scrublet.py
import pandas as pd
import numpy as np
import anndata as ad

# Avoid importing scrublet at top-level elsewhere
import scrublet as scr

h5ad_in = snakemake.input["h5ad"]
tsv_out = snakemake.output["tsv"]
params  = snakemake.params.get("scrublet_params", {}) or {}

adata = ad.read_h5ad(h5ad_in)
# Scrublet expects counts (cells x genes) dense-ish; convert if sparse
X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)

# Allow users to override default params via config
expected_doublet_rate = float(params.get("expected_doublet_rate", 0.06))
min_counts = int(params.get("min_counts", 2))
min_cells = int(params.get("min_cells", 3))

scrub = scr.Scrublet(X, expected_doublet_rate=expected_doublet_rate)
scores, preds = scrub.scrub_doublets(min_counts=min_counts, min_cells=min_cells)

df = pd.DataFrame({
    "barcode": adata.obs_names.astype(str),
    "is_doublet": preds.astype(bool),
    "dbl_score": scores.astype(float),
})
df.to_csv(tsv_out, sep="\t", index=False)
