# lib/api/assemble.py
def assemble_anndata(base_h5ad, artifacts: dict):
    """
    base_h5ad: path to a lean base AnnData (obs/var only, X optional)
    artifacts: dict like {
        "obsm.X_pca": {"path": "results/embeddings/pca/X_pca.npy", "index": "results/embeddings/pca/X_pca.index.csv"},
        "obsp.connectivities": {"path": ".../connectivities.npz"},
        "obs.qc": {"parquet": "results/qc/combined.obs_qc.parquet", "prefix": "raw"},
    }
    """
    import anndata as ad
    adata = ad.read_h5ad(base_h5ad)

    from ..ops.anndata_io import attach_embedding, attach_sparse_obsp, attach_obs_table
    for key, spec in artifacts.items():
        if key.startswith("obsm."):
            attach_embedding(adata, spec["path"], key.split(".",1)[1], spec.get("index"))
        elif key.startswith("obsp."):
            attach_sparse_obsp(adata, spec["path"], key.split(".",1)[1], spec.get("index"))
        elif key == "obs.qc":
            attach_obs_table(adata, spec["parquet"], prefix=spec.get("prefix"))
        # extend as needed

    return adata
