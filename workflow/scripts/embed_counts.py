import anndata as ad
from scpipe.ops.scanpy_ops import pca_from_counts_scanpy
from scpipe.api.embedding import save_embedding_result

adata = ad.read_h5ad(snakemake.input[0])

cfg = dict(
    target_sum = snakemake.params.get("target_sum", 1e4),
    log1p      = snakemake.params.get("log1p", True),
    n_hvgs     = snakemake.params.get("n_hvgs", 2000),
    hvg_flavor = snakemake.params.get("hvg_flavor", "seurat_v3"),
    scale_max  = snakemake.params.get("scale_max", 10.0),
    n_comps    = snakemake.params.get("n_comps", 50),
    zero_center= snakemake.params.get("zero_center", True),
    svd_solver = snakemake.params.get("svd_solver", "arpack"),
    random_state=snakemake.params.get("random_state", 0),
)

res = pca_from_counts_scanpy(adata.X, obs=adata.obs, var=adata.var, **cfg)
save_embedding_result(res, adata.obs_names, snakemake.output.dir)