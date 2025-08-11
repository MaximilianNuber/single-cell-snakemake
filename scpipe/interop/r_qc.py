from __future__ import annotations
from typing import Optional
import numpy as np
import pandas as pd
import anndata as ad

from .r_env import r, get_r_environment
from scpipe.interop.r_env import r
from scpipe.interop.r_env import get_r_environment
from scpipe.interop.r_converters import RConverters as RC

# def run_emptydrops(adata: ad.AnnData, *, lower: int = 100, fdr_cutoff: float = 0.01) -> ad.AnnData:
#     """
#     Run DropletUtils::emptyDrops on counts.
#     Ensures genes x cells orientation, so result index = cell barcodes (matches adata.obs_names).
#     """
#     DropletUtils, = r.lazy_import_r_packages(["DropletUtils"])

#     # genes x cells with correct dimnames
#     r_mat = RC.anndata_to_r_sparse_matrix(adata, orient="genes_by_cells")

#     res = DropletUtils.emptyDrops(r_mat, lower=lower)
#     r_df = r.ro.r["as.data.frame"](res)
#     df = r.r2py(r_df)  # pandas DataFrame, index are cell barcodes (colnames of r_mat)

#     # align back to adata.obs_names
#     df.index = df.index.astype(str)
#     df = df.reindex(adata.obs_names.astype(str))

#     adata.obs["emptydrops_LogProb"] = pd.to_numeric(df.get("LogProb"), errors="coerce")
#     adata.obs["emptydrops_PValue"]  = pd.to_numeric(df.get("PValue"), errors="coerce")
#     adata.obs["emptydrops_FDR"]     = pd.to_numeric(df.get("FDR"), errors="coerce")
#     adata.obs["emptydrops_is_cell"] = (adata.obs["emptydrops_FDR"] < fdr_cutoff).fillna(False).to_numpy()

#     return adata

def run_emptydrops(
    adata: ad.AnnData,
    *,
    lower: int | None = 100,
    fdr_cutoff: float = 0.01,
    retain: int | None = None,
    skip_if_no_ambient: bool = True,
    auto_lower_if_none: bool = True,
    auto_lower_quantile: float = 0.05,
) -> ad.AnnData:
    """
    Run DropletUtils::emptyDrops on counts.
    - If lower is None and auto_lower_if_none, choose it from a counts quantile.
    - If no barcodes fall below 'lower' and skip_if_no_ambient, skip gracefully.
    """
    import numpy as np
    import pandas as pd
    # from scpipe.interop import r_converters as RC, r_env as r  # adjust import to your project

    DropletUtils, = r.lazy_import_r_packages(["DropletUtils"])

    # Determine lower automatically if requested
    if lower is None and auto_lower_if_none:
        n_counts = np.asarray(adata.X.sum(axis=1)).ravel() if hasattr(adata.X, "sum") else adata.X.sum(1)
        # choose a small quantile; clamp to at least 1
        q = int(max(1, np.quantile(n_counts, auto_lower_quantile)))
        lower = q

    # If there are effectively no "ambient" barcodes under lower, optionally skip
    n_counts = np.asarray(adata.X.sum(axis=1)).ravel()
    n_ambient = int((n_counts <= (lower or 0)).sum())
    if skip_if_no_ambient and n_ambient == 0:
        # nothing to test against; return unchanged
        return adata

    # genes x cells with correct dimnames
    r_mat = RC.anndata_to_r_sparse_matrix(adata, orient="genes_by_cells")

    # Call DropletUtils::emptyDrops with lower/retain
    # (retain keeps top barcodes by total count; pass None if not used)
    kwargs = {"lower": int(lower)}
    if retain is not None:
        kwargs["retain"] = int(retain)
    res = DropletUtils.emptyDrops(r_mat, **kwargs)

    r_df = r.ro.r["as.data.frame"](res)
    df = r.r2py(r_df)  # pandas DataFrame, index are cell barcodes (colnames of r_mat)

    df.index = df.index.astype(str)
    df = df.reindex(adata.obs_names.astype(str))

    adata.obs["emptydrops_LogProb"] = pd.to_numeric(df.get("LogProb"), errors="coerce")
    adata.obs["emptydrops_PValue"]  = pd.to_numeric(df.get("PValue"), errors="coerce")
    adata.obs["emptydrops_FDR"]     = pd.to_numeric(df.get("FDR"), errors="coerce")
    adata.obs["emptydrops_is_cell"] = (adata.obs["emptydrops_FDR"] < fdr_cutoff).fillna(False).to_numpy()

    return adata

def run_scdblfinder(adata: ad.AnnData, *, batch_key: Optional[str] = None) -> ad.AnnData:
    """
    scDblFinder on counts. Adds to .obs: is_doublet (bool), dbl_score (float).
    """
    scDblFinder, SingleCellExperiment, SummarizedExperiment = r.lazy_import_r_packages(
        ["scDblFinder", "SingleCellExperiment", "SummarizedExperiment"]
    )
    base = r.ro.baseenv

    # counts genes x cells
    r_mat = RC.anndata_to_r_sparse_matrix(adata, layer=None, orient="genes_by_cells")

    # Build SCE: assays=list(counts=r_mat)
    assays = r.ListVector({"counts": r_mat})
    sce  = SingleCellExperiment.SingleCellExperiment(assays=assays)

    # optional batch
    if batch_key and batch_key in adata.obs:
        sce.colData.rx2["batch"] = r.StrVector(list(map(str, adata.obs[batch_key].astype(str))))

    sce2 = scDblFinder.scDblFinder(sce)
    r_df = r.ro.r["as.data.frame"](SummarizedExperiment.colData(sce2))
    df = r.r2py(r_df)
    df.index = df.index.astype(str)
    df = df.reindex(adata.obs_names.astype(str))

    adata.obs["is_doublet"] = df["scDblFinder.class"].astype(str).eq("doublet").to_numpy()
    adata.obs["dbl_score"]  = pd.to_numeric(df["scDblFinder.score"], errors="coerce").to_numpy()
    return adata

# def run_soupx(adata: ad.AnnData, *, layer_out: str = "soupx_corrected") -> ad.AnnData:
#     """
#     SoupX::autoEstCont + adjCounts. Writes corrected counts into adata.layers[layer_out].
#     Leaves .X as original counts.
#     """
#     SoupX, Matrix, base = r.lazy_import_r_packages(["SoupX", "Matrix", "base"])
#     # Convert counts to genes x cells; prefer sparse dgCMatrix in R for scale
#     # If RC.anndata_to_r_matrix already uses scipy2ri sparse, we're good.
#     r_mat = RC.anndata_to_r_matrix(adata, layer=None)

#     # SoupChannel constructor signatures differ across versions; using tod=, spx= as same here
#     sc = SoupX.SoupChannel(tod=r_mat, spx=r_mat)
#     sc = SoupX.autoEstCont(sc)
#     adj = SoupX.adjCounts(sc)   # genes x cells

#     # back to numpy (cells x genes)
#     adj_np = np.asarray(r.r2py(adj)).T
#     adata.layers[layer_out] = adj_np
#     return adata

def run_soupx_with_raw_filtered(
    adata_filtered: ad.AnnData,
    *,
    raw_adata: Optional[ad.AnnData] = None,
    raw_10x_dir: Optional[str] = None,
    layer_out: str = "soupx_corrected",
    clusters_key: Optional[str] = None,
    orient: Literal["genes_by_cells", "cells_by_genes"] = "genes_by_cells",
) -> ad.AnnData:
    """
    Run SoupX using a raw (droplets) matrix + filtered (cells) matrix.

    Inputs
    ------
    adata_filtered: AnnData
        Cell-called counts (cells x genes). This object will receive the corrected counts in .layers[layer_out].
    raw_adata: Optional[AnnData]
        Unfiltered droplets (cells x genes). If not provided, specify raw_10x_dir.
    raw_10x_dir: Optional[str]
        Path to a 10x 'raw_feature_bc_matrix' directory (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz).
        Ignored if raw_adata is provided.
    layer_out: str
        Name of the target layer to store corrected counts (cells x genes).
    clusters_key: Optional[str]
        If provided and present in adata_filtered.obs, cluster labels are passed to SoupX to improve rho estimation.
    orient: "genes_by_cells" | "cells_by_genes"
        Orientation used when converting to R. SoupX expects genes x cells, so the default is "genes_by_cells".

    Returns
    -------
    adata_filtered (modified in place): corrected counts in .layers[layer_out]
    """

    # --- 0) get raw AnnData ---
    if raw_adata is None and raw_10x_dir is None:
        raise ValueError("Provide either raw_adata or raw_10x_dir for SoupX.")

    if raw_adata is None:
        # lazy import scanpy here to avoid hard dependency in other parts
        import scanpy as sc
        raw_adata = sc.read_10x_mtx(raw_10x_dir, var_names="gene_symbols", cache=False)
        raw_adata.var_names_make_unique()

    # --- 1) subset raw to genes present in filtered, and keep all droplets ---
    # Make sure gene names align across both objects
    common_genes = adata_filtered.var_names.intersection(raw_adata.var_names)
    if len(common_genes) == 0:
        raise ValueError("No overlapping genes between filtered and raw matrices for SoupX.")

    raw_use = raw_adata[:, common_genes].copy()
    filt_use = adata_filtered[:, common_genes].copy()

    # --- 2) build R sparse matrices with correct orientation (SoupX expects genes x cells) ---
    r_raw = RC.anndata_to_r_sparse_matrix(raw_use, orient=orient)   # genes x droplets
    r_flt = RC.anndata_to_r_sparse_matrix(filt_use, orient=orient)  # genes x cells

    # --- 3) create SoupChannel (different SoupX versions use tod+toc or tod+spx) ---
    SoupX, = r.lazy_import_r_packages(["SoupX"])

    sc_obj = None
    # Try tod/toc (most common)
    try:
        sc_obj = SoupX.SoupChannel(tod=r_raw, toc=r_flt)
    except Exception:
        # Fallback: some versions accept 'spx' instead of 'toc'
        sc_obj = SoupX.SoupChannel(tod=r_raw, spx=r_flt)

    # --- 4) (optional) provide clusters to SoupX (improves rho) ---
    if clusters_key and clusters_key in adata_filtered.obs.columns:
        # Ensure vector is ordered like filtered colnames; our r_flt has colnames set to adata_filtered.obs_names
        cl = pd.Series(adata_filtered.obs[clusters_key].astype(str), index=adata_filtered.obs_names.astype(str))
        # Align to R colnames on sc_obj
        r_colnames = list(r.ro.baseenv["colnames"](r_flt))
        cl = cl.reindex(r_colnames)
        sc_obj = SoupX.setClusters(sc_obj, r.StrVector(list(map(str, cl.fillna("NA")))))

    # --- 5) estimate contamination and adjust counts ---
    sc_obj = SoupX.autoEstCont(sc_obj)
    # Different versions: adjustCounts vs adjCounts
    try:
        adj = SoupX.adjustCounts(sc_obj)
    except Exception:
        adj = SoupX.adjCounts(sc_obj)

    # adj is genes x cells; convert back to numpy (cells x genes) with the same cell order as filtered
    adj_np = np.asarray(r.r2py(adj)).T  # cells x genes

    # --- 6) Assign corrected counts back to adata_filtered in the same cell order ---
    # Ensure columns (cells) from adj match obs_names (they should, since we set colnames earlier)
    r_cols = np.asarray(r.ro.baseenv["colnames"](adj)).astype(str)
    if not np.array_equal(r_cols, adata_filtered.obs_names.astype(str).to_numpy()):
        # align safely in case R reorders
        pos = pd.Index(r_cols).get_indexer(adata_filtered.obs_names.astype(str))
        if np.any(pos < 0):
            missing = set(adata_filtered.obs_names.astype(str)) - set(r_cols)
            raise RuntimeError(f"SoupX returned matrix missing some filtered cells (e.g., {list(missing)[:5]})")
        adj_np = adj_np[pos, :]

    # Ensure gene order
    r_rows = np.asarray(r.ro.baseenv["rownames"](adj)).astype(str)
    if not np.array_equal(r_rows, common_genes.astype(str).to_numpy()):
        # Align to common_genes order
        posg = pd.Index(r_rows).get_indexer(common_genes.astype(str))
        adj_np = adj_np[:, posg]

    # Attach layer (on the full filtered AnnData; genes not in common_genes will be absent)
    # To keep shape consistent with adata_filtered (all genes), we can create a zeros array and fill the common block.
    out = np.zeros((adata_filtered.n_obs, adata_filtered.n_vars), dtype=adj_np.dtype)
    gi = adata_filtered.var_names.get_indexer(common_genes)
    out[:, gi] = adj_np
    adata_filtered.layers[layer_out] = out

    # store a bit of metadata
    adata_filtered.uns.setdefault("soupx", {})["layer"] = layer_out
    return adata_filtered