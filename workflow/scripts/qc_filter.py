# workflow/scripts/qc_filter.py
import pandas as pd
import anndata as ad
from scpipe.ops.qc_scanpy import QCFilterConfig, add_mad_outlier_calls, apply_qc_filters

h5ad_in = snakemake.input["h5ad"]
dbl_tsv = snakemake.input["dbl"]
h5ad_out = snakemake.output["h5ad"]
qc_cfg_dict = snakemake.params["qc_cfg"]
percent_top = int(snakemake.params["percent_top"])

adata = ad.read_h5ad(h5ad_in)

# Merge doublets if a TSV was produced (else it's just header)
if pd.read_csv(dbl_tsv, sep="\t").shape[0] > 0:
    df = pd.read_csv(dbl_tsv, sep="\t")
    df.set_index("barcode", inplace=True)
    df = df.reindex(adata.obs_names.astype(str))
    if "is_doublet" in df:
        adata.obs["is_doublet"] = df["is_doublet"].fillna(False).astype(bool).to_numpy()
    if "dbl_score" in df:
        adata.obs["dbl_score"] = pd.to_numeric(df["dbl_score"], errors="coerce").to_numpy()

# Create filter config
cfg = QCFilterConfig(
    mad_total_counts=qc_cfg_dict.get("nmads_total_counts", 5),
    mad_n_genes=qc_cfg_dict.get("nmads_n_genes", 5),
    mad_pct_top=qc_cfg_dict.get("nmads_pct_top", 5),
    mad_pct_mito=qc_cfg_dict.get("nmads_pct_mito", 3),
    max_pct_mito=qc_cfg_dict.get("max_pct_mito", 8),
    min_cells_per_gene=qc_cfg_dict.get("min_cells_per_gene", 3),
    min_genes_per_cell=qc_cfg_dict.get("min_genes_per_cell", None),
    keep_doublets=qc_cfg_dict.get("keep_doublets", False),
)

add_mad_outlier_calls(adata, cfg, percent_top_bucket=percent_top)
adata_filt = apply_qc_filters(adata, cfg)
adata_filt.write_h5ad(h5ad_out, compression="gzip")
