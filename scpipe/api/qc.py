# workflow/lib/api/qc.py
from __future__ import annotations
from pathlib import Path
from typing import Optional, Literal, Mapping
import anndata as ad
import scanpy as sc

from ..ops.qc_scanpy import (
    QCMetricsConfig, QCFilterConfig,
    compute_qc_metrics, add_mad_outlier_calls, apply_qc_filters,
)
from ..interop.r_qc import run_emptydrops, run_scdblfinder, run_soupx_with_raw_filtered
from .qc_configs import EmptyDropsConfig, ScDblFinderConfig, SoupXConfig


def read_cellranger(matrix: str | Path, var_names: Literal["gene_symbols","gene_ids"]="gene_symbols") -> ad.AnnData:
    p = Path(matrix)
    if p.is_dir():
        adata = sc.read_10x_mtx(p, var_names=var_names, cache=False)
    else:
        adata = sc.read_10x_h5(p)
        if var_names == "gene_symbols" and "gene_symbols" in adata.var:
            adata.var_names = adata.var["gene_symbols"].astype(str)
    adata.var_names_make_unique()
    return adata

def _maybe_sibling_raw_dir(filtered_dir: Path) -> Optional[Path]:
    if filtered_dir.is_dir() and filtered_dir.name == "filtered_feature_bc_matrix":
        candidate = filtered_dir.parent / "raw_feature_bc_matrix"
        return candidate if candidate.is_dir() else None
    return None


####

def _coerce_emptydrops(x: Optional[Union[EmptyDropsConfig, Mapping]]) -> Optional[EmptyDropsConfig]:
    if x is None or isinstance(x, EmptyDropsConfig):
        return x
    if isinstance(x, Mapping):
        return EmptyDropsConfig(
            lower=int(x.get("lower", 100)),
            fdr_cutoff=float(x.get("fdr_cutoff", x.get("fdr", 0.01))),
            skip_if_no_ambient=bool(x.get("skip_if_no_ambient", True)),
            auto_lower_if_none=bool(x.get("auto_lower_if_none", True)),
            auto_lower_quantile=float(x.get("auto_lower_quantile", 0.05)),
        )
    raise TypeError(f"emptydrops must be EmptyDropsConfig or dict, got {type(x).__name__}")

def _coerce_soupx(x: Optional[Union[SoupXConfig, Mapping]]) -> Optional[SoupXConfig]:
    if x is None or isinstance(x, SoupXConfig):
        return x
    if isinstance(x, Mapping):
        return SoupXConfig(
            raw_10x_dir=x.get("raw_10x_dir"),
            layer_out=x.get("layer_out", "soupx"),
            clusters_key=x.get("clusters_key", None),
        )
    raise TypeError(f"soupx must be SoupXConfig or dict, got {type(x).__name__}")

def _coerce_scdbl(x: Optional[Union[ScDblFinderConfig, Mapping]]) -> Optional[ScDblFinderConfig]:
    if x is None or isinstance(x, ScDblFinderConfig):
        return x
    if isinstance(x, Mapping):
        return ScDblFinderConfig(batch_key=x.get("batch_key", None))
    raise TypeError(f"scdblfinder must be ScDblFinderConfig or dict, got {type(x).__name__}")

####
def qc_compute_metrics(
    adata: ad.AnnData,
    *,
    metrics_cfg: QCMetricsConfig = QCMetricsConfig(),
    emptydrops: Optional[EmptyDropsConfig] = None,
    scdblfinder: Optional[ScDblFinderConfig] = None,
    soupx: Optional[SoupXConfig] = None,
    filtered_matrix_path: Optional[Path] = None,  # for raw sibling discovery
) -> ad.AnnData:
    """Compute Scanpy QC metrics and optional R-based tools."""

    emptydrops  = _coerce_emptydrops(emptydrops)
    scdblfinder = _coerce_scdbl(scdblfinder)
    soupx       = _coerce_soupx(soupx)


    compute_qc_metrics(adata, metrics_cfg)

    if emptydrops is not None:
        adata = run_emptydrops(
            adata,
            lower=emptydrops.lower,
            fdr_cutoff=emptydrops.fdr_cutoff,
            skip_if_no_ambient=emptydrops.skip_if_no_ambient,
            auto_lower_if_none=emptydrops.auto_lower_if_none,
            auto_lower_quantile=emptydrops.auto_lower_quantile,
        )

    if scdblfinder is not None:
        adata = run_scdblfinder(
            adata,
            batch_key=scdblfinder.batch_key,
        )

    if soupx is not None:
        raw_dir = soupx.raw_10x_dir
        if raw_dir is None and filtered_matrix_path is not None:
            sib = _maybe_sibling_raw_dir(Path(filtered_matrix_path))
            if sib is not None:
                raw_dir = str(sib)
        if raw_dir is None:
            raise ValueError("SoupXConfig requires raw_10x_dir (or pass filtered_matrix_path so we can infer sibling raw).")

        adata = run_soupx_with_raw_filtered(
            adata_filtered=adata,
            raw_10x_dir=raw_dir,
            layer_out=soupx.layer_out,
            clusters_key=soupx.clusters_key,
        )

    return adata

def qc_filter(
    adata: ad.AnnData,
    *,
    filter_cfg: QCFilterConfig = QCFilterConfig(),
    percent_top_bucket: int = 20,
) -> ad.AnnData:
    add_mad_outlier_calls(adata, filter_cfg, percent_top_bucket=percent_top_bucket)
    return apply_qc_filters(adata, filter_cfg)

def qc_single_sample_from_cellranger(
    *,
    matrix: str | Path,
    sample_name: Optional[str] = None,
    var_names: Literal["gene_symbols","gene_ids"]="gene_symbols",
    metrics_cfg: QCMetricsConfig = QCMetricsConfig(),
    filter_cfg: QCFilterConfig = QCFilterConfig(),
    percent_top_bucket: int = 20,
    emptydrops: Optional[EmptyDropsConfig] = None,
    scdblfinder: Optional[ScDblFinderConfig] = None,
    soupx: Optional[SoupXConfig] = None,
    out_h5ad_metrics: Optional[str | Path] = None,
    out_h5ad_filtered: Optional[str | Path] = None,
) -> tuple[ad.AnnData, ad.AnnData]:
    # Load filtered (or raw) matrix
    adata = read_cellranger(matrix, var_names=var_names)
    if sample_name:
        adata.obs["sample"] = str(sample_name)

    # Compute metrics + optional tools
    filtered_matrix_path = Path(matrix) if Path(matrix).is_dir() else None
    qc_compute_metrics(
        adata,
        metrics_cfg=metrics_cfg,
        emptydrops=emptydrops,
        scdblfinder=scdblfinder,
        soupx=soupx,
        filtered_matrix_path=filtered_matrix_path,
    )

    if out_h5ad_metrics:
        Path(out_h5ad_metrics).parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(out_h5ad_metrics, compression="gzip")

    # Filter
    adata_filt = qc_filter(adata, filter_cfg=filter_cfg, percent_top_bucket=percent_top_bucket)

    if out_h5ad_filtered:
        Path(out_h5ad_filtered).parent.mkdir(parents=True, exist_ok=True)
        adata_filt.write_h5ad(out_h5ad_filtered, compression="gzip")

    return adata, adata_filt



# from __future__ import annotations
# from pathlib import Path
# from typing import Optional, Literal
# import anndata as ad
# import scanpy as sc

# from ..ops.qc_scanpy import (
#     QCMetricsConfig, QCFilterConfig,
#     compute_qc_metrics, add_mad_outlier_calls, apply_qc_filters,
# )
# from ..interop.r_qc import run_emptydrops, run_scdblfinder, run_soupx_with_raw_filtered
# from .qc_configs import EmptyDropsConfig, ScDblFinderConfig, SoupXConfig

# def read_cellranger(matrix: str | Path, var_names: Literal["gene_symbols","gene_ids"]="gene_symbols") -> ad.AnnData:
#     p = Path(matrix)
#     if p.is_dir():
#         adata = sc.read_10x_mtx(p, var_names=var_names, cache=False)
#     else:
#         adata = sc.read_10x_h5(p)
#         if var_names == "gene_symbols" and "gene_symbols" in adata.var:
#             adata.var_names = adata.var["gene_symbols"].astype(str)
#     adata.var_names_make_unique()
#     return adata

# def qc_compute_metrics(
#     adata: ad.AnnData,
#     *,
#     metrics_cfg: QCMetricsConfig = QCMetricsConfig(),
#     run_emptydrops_fdr: Optional[float] = None,
#     run_soupx_layer_out: Optional[str] = None,
#     run_scdblfinder_batch: Optional[str] = None,
# ) -> ad.AnnData:
#     compute_qc_metrics(adata, metrics_cfg)
#     if run_emptydrops_fdr is not None:
#         adata = run_emptydrops(adata, fdr_cutoff=float(run_emptydrops_fdr))
#     if run_scdblfinder_batch is not None:
#         adata = run_scdblfinder(adata, batch_key=run_scdblfinder_batch)
#     if run_soupx_layer_out is not None:
#         adata = run_soupx(adata, layer_out=run_soupx_layer_out)
#     return adata

# def qc_filter(
#     adata: ad.AnnData,
#     *,
#     filter_cfg: QCFilterConfig = QCFilterConfig(),
#     percent_top_bucket: int = 20,
# ) -> ad.AnnData:
#     add_mad_outlier_calls(adata, filter_cfg, percent_top_bucket=percent_top_bucket)
#     return apply_qc_filters(adata, filter_cfg)

# def qc_single_sample_from_cellranger(
#     *,
#     matrix: str | Path,
#     sample_name: Optional[str] = None,
#     var_names: Literal["gene_symbols","gene_ids"]="gene_symbols",
#     metrics_cfg: QCMetricsConfig = QCMetricsConfig(),
#     filter_cfg: QCFilterConfig = QCFilterConfig(),
#     percent_top_bucket: int = 20,
#     run_emptydrops_fdr: Optional[float] = None,
#     run_soupx_layer_out: Optional[str] = None,
#     run_scdblfinder_batch: Optional[str] = None,
#     out_h5ad_metrics: Optional[str | Path] = None,
#     out_h5ad_filtered: Optional[str | Path] = None,
# ):
#     adata = read_cellranger(matrix, var_names=var_names)
#     if sample_name: adata.obs["sample"] = str(sample_name)
#     qc_compute_metrics(
#         adata,
#         metrics_cfg=metrics_cfg,
#         run_emptydrops_fdr=run_emptydrops_fdr,
#         run_soupx_layer_out=run_soupx_layer_out,
#         run_scdblfinder_batch=run_scdblfinder_batch,
#     )
#     if out_h5ad_metrics:
#         Path(out_h5ad_metrics).parent.mkdir(parents=True, exist_ok=True)
#         adata.write_h5ad(out_h5ad_metrics, compression="gzip")
#     adata_filt = qc_filter(adata, filter_cfg=filter_cfg, percent_top_bucket=percent_top_bucket)
#     if out_h5ad_filtered:
#         Path(out_h5ad_filtered).parent.mkdir(parents=True, exist_ok=True)
#         adata_filt.write_h5ad(out_h5ad_filtered, compression="gzip")
#     return adata, adata_filt
