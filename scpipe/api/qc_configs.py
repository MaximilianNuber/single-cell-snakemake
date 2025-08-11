from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

@dataclass(frozen=True)
class EmptyDropsConfig:
    lower: int = 100
    fdr_cutoff: float = 0.01
    skip_if_no_ambient: bool = True
    auto_lower_if_none: bool = False
    auto_lower_quantile: float = 0.01  # 1%

@dataclass(frozen=True)
class ScDblFinderConfig:
    batch_key: Optional[str] = None

@dataclass(frozen=True)
class SoupXConfig:
    # You must provide one of these:
    raw_10x_dir: Optional[str] = None
    # or, if you already loaded it:
    # raw_adata: Optional[anndata.AnnData] = None   # (we’ll pass this at call site if needed)
    layer_out: str = "soupx_corrected"
    clusters_key: Optional[str] = None