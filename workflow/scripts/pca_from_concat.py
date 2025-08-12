from __future__ import annotations
from pathlib import Path
import anndata as ad
from scpipe.ops.pca_ops import pca_from_counts_scanpy
from scpipe.api.embedding import save_embedding_result

adata = ad.read_h5ad(snakemake.input[0])
cfg = dict(snakemake.params.get("pca_cfg", {}))

res = pca_from_counts_scanpy(adata, **cfg)

# Support both styles: named file outputs or a directory output
if isinstance(snakemake.output, dict) and "x" in snakemake.output:
    out_dir = Path(snakemake.output["x"]).parent
else:
    # assume first (or only) output is the directory
    out_dir = Path(next(iter(snakemake.output)))

out_dir.mkdir(parents=True, exist_ok=True)
save_embedding_result(res, adata.obs_names, out_dir, save_loadings_=False)