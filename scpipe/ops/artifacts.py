# lib/ops/artifacts.py
from __future__ import annotations
import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Sequence

def save_embedding(
    X: np.ndarray,
    obs_names: Sequence[str],
    out_dir: str | Path,
    key: str = "X_pca",
) -> dict:
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    p_npy = out / f"{key}.npy"
    p_idx = out / f"{key}.index.csv"
    np.save(p_npy, X)
    pd.Series(obs_names, name="obs_names").to_csv(p_idx, index=False)
    return {"path": str(p_npy), "index": str(p_idx)}

def save_loadings(
    PCs: np.ndarray,
    var_names: Sequence[str],
    out_dir: str | Path,
    key: str = "PCs",
) -> dict:
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    p_npy = out / f"{key}.npy"
    p_idx = out / f"{key}.genes.csv"
    np.save(p_npy, PCs)
    pd.Series(var_names, name="var_names").to_csv(p_idx, index=False)
    return {"path": str(p_npy), "genes": str(p_idx)}

def _to_jsonable(x):
    # numpy scalars
    if isinstance(x, np.generic):
        return x.item()
    # numpy arrays
    if isinstance(x, np.ndarray):
        return x.tolist()
    # pandas
    if isinstance(x, pd.Series):
        return x.astype(object).where(pd.notna(x), None).tolist()
    if isinstance(x, pd.DataFrame):
        return x.to_dict(orient="list")
    return x  # let json handle dict/list/str/float/None

def save_json(obj: dict, out_path: str | Path) -> str:
    out = Path(out_path); out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump(obj, f, indent=2, default=_to_jsonable)
    return str(out)

def save_sparse_csr(csr, out_path):
    from scipy.sparse import save_npz
    from pathlib import Path
    out = Path(out_path); out.parent.mkdir(parents=True, exist_ok=True)
    save_npz(out, csr)
    return str(out)

def save_labels(labels, obs_names, out_csv, key="leiden"):
    import pandas as pd
    from pathlib import Path
    out = Path(out_csv); out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame({key: labels}, index=pd.Index(obs_names, name="obs_names"))
    df.to_csv(out)
    return str(out)