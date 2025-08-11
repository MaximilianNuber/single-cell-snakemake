# workflow/scripts/qc_compute_metrics.py
import json
import anndata as ad
from scpipe.api.qc import read_cellranger, qc_compute_metrics, QCMetricsConfig
from scpipe.api.qc_configs import EmptyDropsConfig, SoupXConfig, ScDblFinderConfig  # optional usage

matrix = snakemake.input["matrix"]
out    = snakemake.output["h5ad"]
percent_top = snakemake.params["percent_top"]
# After (safe)
# emptydrops_cfg = snakemake.config.get("emptydrops", {})
# soupx_cfg     = snakemake.config.get("soupx", {})
# scdbl_cfg     = snakemake.config.get("scdbl", {})  # NOTE: we won't run scDblFinder here; separate rule below

# Read from the correct place and default to None (not {})
tools = snakemake.config.get("tools", {})
emptydrops_cfg = tools.get("emptydrops")    # None if null/missing
soupx_cfg      = tools.get("soupx")         # None if null/missing

# pick up scDblFinder only if selected
# dbl = tools.get("doublets", {})
# scdbl_cfg = None
# if dbl.get("method") == "scdblfinder":
#     p = dbl.get("params", {}) or {}
#     scdbl_cfg = ScDblFinderConfig(batch_key=p.get("batch_key"))

adata = read_cellranger(matrix, var_names="gene_symbols")
cfg = QCMetricsConfig(percent_top=tuple(percent_top))

# Only Scanpy metrics here, keep R-based steps to their own rules if you prefer strict env isolation.
adata = qc_compute_metrics(
    adata,
    metrics_cfg=cfg,
    emptydrops=emptydrops_cfg,    # set to EmptyDropsConfig(...) if you want to run it here with r_qc.yaml
    soupx=soupx_cfg,         # same note
    scdblfinder=None,   # scDblFinder runs in its own rule
    filtered_matrix_path=None,
)

adata.write_h5ad(out, compression="gzip")
