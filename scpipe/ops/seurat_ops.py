# lib/ops/seurat_ops.py
from ..types import Result

def seurat_lognormalize(counts_df, scale_factor: float = 1e4) -> Result:
    """
    counts_df: pandas DataFrame (cells x genes)
    Returns: normalized matrix as outputs["layer.norm"]
    """
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()

    ro.r("suppressPackageStartupMessages(library(Seurat))")
    r_counts = pandas2ri.py2rpy(counts_df)
    ro.globalenv["m"] = r_counts
    ro.r(f"obj <- CreateSeuratObject(counts = t(m))")  # Seurat expects genes x cells
    ro.r(f"obj <- NormalizeData(obj, normalization.method = 'LogNormalize', scale.factor = {scale_factor})")
    norm = ro.r("t(as.matrix(obj[['RNA']]@data))")     # back to cells x genes
    norm_df = pandas2ri.rpy2py(norm)
    return Result(outputs={"layer.norm": norm_df.values}, state={"scale_factor": scale_factor})