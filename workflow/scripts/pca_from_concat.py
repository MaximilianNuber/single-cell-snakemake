# workflow/scripts/pca_from_concat.py
import anndata as ad
from scpipe.ops.scanpy_pca import pca_from_counts_scanpy
from scpipe.api.embedding import save_embedding_result

adata = ad.read_h5ad(snakemake.input[0])
cfg = dict(snakemake.params.get("pca_cfg", {}))

res = pca_from_counts_scanpy(adata, **cfg)  # accepts AnnData directly
save_embedding_result(res, adata.obs_names, snakemake.output[0], save_loadings_=False)
