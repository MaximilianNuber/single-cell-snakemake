from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence, Optional
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.stats import median_abs_deviation

@dataclass(frozen=True)
class QCMetricsConfig:
    mito_prefixes: Sequence[str] = ("MT-", "mt-")
    ribo_prefixes: Sequence[str] = ("RPS", "RPL", "Mrps", "Mrpl")
    hb_regex: str = r"^HB[^(P)]"
    percent_top: Sequence[int] = (20,)
    log1p: bool = True

@dataclass(frozen=True)
class QCFilterConfig:
    mad_total_counts: float = 5.0
    mad_n_genes: float = 5.0
    mad_pct_top: float = 5.0
    mad_pct_mito: float = 3.0
    max_pct_mito: Optional[float] = 8.0
    min_cells_per_gene: int = 3
    min_genes_per_cell: Optional[int] = None
    keep_doublets: bool = False

def _flag_genes(adata, cfg):
    vn = pd.Index(adata.var_names.astype(str))

    mito = np.zeros(adata.n_vars, dtype=bool)
    for pre in cfg.mito_prefixes:
        mito |= np.asarray(vn.str.startswith(pre, na=False))

    ribo = np.zeros(adata.n_vars, dtype=bool)
    for pre in cfg.ribo_prefixes:
        ribo |= np.asarray(vn.str.startswith(pre, na=False))

    hb = np.asarray(vn.str.contains(cfg.hb_regex, regex=True, na=False))

    adata.var["mt"] = mito
    adata.var["ribo"] = ribo
    adata.var["hb"] = hb

def _mad_outlier(s: pd.Series | np.ndarray, nmads: float) -> np.ndarray:
    x = pd.Series(s, dtype=float)
    med = np.nanmedian(x)
    mad = median_abs_deviation(x, nan_policy="omit", scale="normal")
    if mad == 0 or np.isnan(mad): return np.zeros(len(x), bool)
    lo, hi = med - nmads * mad, med + nmads * mad
    return (x < lo) | (x > hi)

def compute_qc_metrics(adata: ad.AnnData, cfg: QCMetricsConfig = QCMetricsConfig()) -> ad.AnnData:
    _flag_genes(adata, cfg)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt","ribo","hb"],
        percent_top=tuple(cfg.percent_top),
        log1p=cfg.log1p,
        inplace=True,
    )
    return adata

def add_mad_outlier_calls(adata: ad.AnnData, cfg: QCFilterConfig = QCFilterConfig(), percent_top_bucket: int = 20):
    key = f"pct_counts_in_top_{percent_top_bucket}_genes"
    if key not in adata.obs: raise KeyError(f"Run compute_qc_metrics with percent_top including {percent_top_bucket}.")
    adata.obs["outlier"] = (
        _mad_outlier(adata.obs["log1p_total_counts"], cfg.mad_total_counts)
        | _mad_outlier(adata.obs["log1p_n_genes_by_counts"], cfg.mad_n_genes)
        | _mad_outlier(adata.obs[key], cfg.mad_pct_top)
    )
    mt_mad = _mad_outlier(adata.obs["pct_counts_mt"], cfg.mad_pct_mito)
    mt_cap = (adata.obs["pct_counts_mt"] > cfg.max_pct_mito).to_numpy() if cfg.max_pct_mito is not None else False
    adata.obs["mt_outlier"] = mt_mad | mt_cap
    adata.obs["qc_fail"] = adata.obs["outlier"] | adata.obs["mt_outlier"]
    return adata

def apply_qc_filters(adata: ad.AnnData, cfg: QCFilterConfig = QCFilterConfig()) -> ad.AnnData:
    keep = ~adata.obs.get("qc_fail", pd.Series(False, index=adata.obs_names)).astype(bool).to_numpy()
    if not cfg.keep_doublets and "is_doublet" in adata.obs: keep &= ~adata.obs["is_doublet"].astype(bool).to_numpy()
    out = adata[keep].copy()
    if cfg.min_genes_per_cell is not None: sc.pp.filter_cells(out, min_genes=cfg.min_genes_per_cell)
    if cfg.min_cells_per_gene and cfg.min_cells_per_gene > 0: sc.pp.filter_genes(out, min_cells=cfg.min_cells_per_gene)
    return out
