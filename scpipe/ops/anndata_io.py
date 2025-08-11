# lib/ops/anndata_io.py

from __future__ import annotations
def write_layer(adata, arr, name: str): adata.layers[name] = arr
def write_obs(adata, key: str, values): adata.obs[key] = values
def write_obsm(adata, key: str, arr):   adata.obsm[key] = arr


# save
def save_embedding(X_pca, obs_names, path_npy, path_index):
    np.save(path_npy, X_pca)
    pd.Series(obs_names, name="obs_names").to_csv(path_index, index=False)

def save_sparse_csr(csr, path_npz, path_index):
    from scipy.sparse import save_npz
    import pandas as pd
    save_npz(path_npz, csr)
    pd.Series(csr.shape, index=["n_rows","n_cols"]).to_csv(path_index)

def save_obs_table(df, path_parquet):
    df.to_parquet(path_parquet)

# load & attach
def attach_embedding(adata, path_npy, key, path_index=None):
    import numpy as np, pandas as pd
    X = np.load(path_npy)
    if path_index:
        idx = pd.read_csv(path_index)["obs_names"].astype(str).to_numpy()
        X = X[_align_index(idx, adata.obs_names)]
    adata.obsm[key] = X
    return adata

def attach_sparse_obsp(adata, path_npz, key, path_index=None):
    from scipy.sparse import load_npz
    M = load_npz(path_npz)
    adata.obsp[key] = M
    return adata

def attach_obs_table(adata, path_parquet, prefix=None):
    import pandas as pd
    df = pd.read_parquet(path_parquet)
    df = df.reindex(adata.obs_names)
    if prefix: df = df.add_prefix(prefix + "_")
    for c in df.columns: adata.obs[c] = df[c].values
    return adata

def _align_index(source_names, target_names):
    import numpy as np
    pos = {name: i for i, name in enumerate(source_names)}
    return np.array([pos[n] for n in target_names], dtype=int)
