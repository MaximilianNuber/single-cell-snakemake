import os
from typing import Any, Optional, Sequence, Literal

import pandas as pd
import numpy as np
import anndata 
import scipy.sparse as sp

from .r_env import get_r_environment, r

# A type hint for an AnnData object
AnnData = anndata.AnnData

class RConverters:
    """A namespace for custom rpy2 conversion methods."""

    @staticmethod
    def pandas_to_r_matrix(df: pd.DataFrame) -> Any:
        """Converts a pandas DataFrame to an R matrix."""
        r = get_r_environment()
        with r.localconverter(r.default_converter + r.numpy2ri_converter):
            r_mat = r.get_conversion().py2rpy(df.values)
        # dimnames
        r_mat = r.ro.baseenv["rownames<-"](r_mat, r.StrVector(df.index.astype(str).to_numpy()))
        r_mat = r.ro.baseenv["colnames<-"](r_mat, r.StrVector(df.columns.astype(str).to_numpy()))
        return r_mat

    @staticmethod
    def pandas_to_r_df(df: pd.DataFrame) -> Any:
        """Converts a pandas DataFrame to an R data.frame."""
        r = get_r_environment()
        with r.localconverter(r.default_converter + r.pandas2ri_converter):
            return r.get_conversion().py2rpy(df)

    @staticmethod
    def numpy_to_r_matrix(arr: np.ndarray, rownames: Optional[Sequence[str]] = None,
                       colnames: Optional[Sequence[str]] = None) -> Any:
        """Converts a numpy array to an R matrix."""
        r = get_r_environment()
        with r.localconverter(r.default_converter + r.numpy2ri_converter):
            r_mat = r.get_conversion().py2rpy(arr)
        if rownames is not None:
            r_mat = r.ro.baseenv["rownames<-"](r_mat, r.StrVector(list(map(str, rownames))))
        if colnames is not None:
            r_mat = r.ro.baseenv["colnames<-"](r_mat, r.StrVector(list(map(str, colnames))))
        return r_mat

    @staticmethod
    def anndata_to_r_matrix(adata: AnnData, layer: Optional[str] = None) -> Any:
        """
        Converts AnnData count matrix to an R matrix (genes x samples).
        """
        r = get_r_environment()
        import scipy.sparse as sp
        mat = adata.X if layer is None else adata.layers[layer]
        if sp.issparse(mat):
            with r.localconverter(r.default_converter + r.scipy2ri.converter):
                r_mat = r.get_conversion().py2rpy(mat)
        else:
            with r.localconverter(r.default_converter + r.numpy2ri_converter):
                r_mat = r.get_conversion().py2rpy(mat)
        # Set rownames and colnames
        if adata.var_names is not None:
            r_mat = r.ro.baseenv["colnames<-"](r_mat, r.StrVector(adata.var_names.astype(str).to_numpy()))
        if adata.obs_names is not None:
            r_mat = r.ro.baseenv["rownames<-"](r_mat, r.StrVector(adata.obs_names.astype(str).to_numpy()))
        
        return r_mat
    
    @staticmethod
    def rmatrix_to_pandas(r_mat: Any) -> pd.DataFrame:
        """Converts an R matrix to a pandas DataFrame."""
        r = get_r_environment()
        with r.localconverter(r.default_converter + r.numpy2ri_converter):
            pymat = r.get_conversion().rpy2py(r_mat)
        df = pd.DataFrame(pymat)
        # set rownames and colnames

        colnames = np.asarray(
            r.ro.baseenv["colnames"](r_mat), dtype=str
        )
        rownames = np.asarray(
            r.ro.baseenv["rownames"](r_mat), dtype=str
        )
        if colnames is not None:
            df.columns = colnames
        if rownames is not None:
            df.index = rownames

        return df

    @staticmethod
    def rmatrix_to_numpy(r_mat: Any) -> np.ndarray:
        """Converts an R matrix to a numpy array."""
        r = get_r_environment()
        with r.localconverter(r.default_converter + r.numpy2ri_converter):
            return r.get_conversion().rpy2py(r_mat)

    @staticmethod
    def rmatrix_to_anndata(r_mat: Any) -> AnnData:
        """
        Converts an R matrix (genes x samples) to an AnnData object (samples x genes).
        """
        # Convert R matrix to numpy, get dimnames
        r = get_r_environment()
        df = _RConverters.rmatrix_to_numpy(r_mat)
        rcolnames = r.ro.baseenv["colnames"](r_mat)
        rrownames = r.ro.baseenv["rownames"](r_mat)

        print(np.asarray(rcolnames))
        print(isinstance(np.asarray(rcolnames), type(r.ro.NULL)))

        colnames = None if isinstance((rcolnames), type(r.ro.NULL)) else np.asarray(rcolnames)
        rownames = None if isinstance((rrownames), type(r.ro.NULL)) else np.asarray(rrownames)

        adata = AnnData(df.T)
        print(colnames)
        print(rownames)
        if colnames is not None:
            adata.var_names = rownames
        if rownames is not None:
            adata.obs_names = colnames

        return adata
    
    # ---------------- NEW: sparse -> R sparse ----------------
    @staticmethod
    def sparse_matrix_to_r_sparse_matrix(mat: sp.spmatrix) -> Any:
        """
        Convert a SciPy sparse matrix to an R Matrix::dgCMatrix (column-compressed).
        Never densifies. Uses anndata2ri/scipy2ri bridge from your r_env.
        """
        r = get_r_environment()
        if not sp.issparse(mat):
            # be strict: caller should pass a sparse matrix;
            # if you want to be permissive, uncomment next line
            # mat = sp.csc_matrix(mat)
            raise TypeError("Expected a SciPy sparse matrix (csr/csc/coo).")

        # R's dgCMatrix is column-compressed; prefer csc on the Python side
        m_csc = mat.tocsc(copy=False)
        with r.localconverter(r.default_converter + r.scipy2ri.converter):
            r_mat = r.get_conversion().py2rpy(m_csc)
        return r_mat

    # ---------------- NEW: AnnData -> R sparse with names ----------------
    @staticmethod
    def anndata_to_r_sparse_matrix(
        adata: AnnData,
        layer: Optional[str] = None,
        *,
        orient: Literal["genes_by_cells", "cells_by_genes"] = "genes_by_cells",
    ) -> Any:
        """
        Convert AnnData counts to an R Matrix::dgCMatrix with correct orientation and dimnames.

        orient="genes_by_cells"  -> rows = genes (var_names), cols = cells (obs_names) [default]
        orient="cells_by_genes"  -> rows = cells (obs_names), cols = genes (var_names)
        """
        r = get_r_environment()
        X = adata.X if layer is None else adata.layers[layer]

        if not sp.issparse(X):
            # keep memory-friendly; if you prefer strictness, raise instead
            X = sp.csr_matrix(X)

        # Shape in AnnData is (cells x genes). Choose orientation:
        if orient == "genes_by_cells":
            pymat = X.T   # genes x cells
            rownames = adata.var_names.astype(str).to_numpy()
            colnames = adata.obs_names.astype(str).to_numpy()
        elif orient == "cells_by_genes":
            pymat = X     # cells x genes
            rownames = adata.obs_names.astype(str).to_numpy()
            colnames = adata.var_names.astype(str).to_numpy()
        else:
            raise ValueError("orient must be 'genes_by_cells' or 'cells_by_genes'.")

        # Convert sparse -> R dgCMatrix
        r_mat = RConverters.sparse_matrix_to_r_sparse_matrix(pymat)

        # Assign dimnames
        r_mat = r.ro.baseenv["rownames<-"](r_mat, r.StrVector(list(map(str, rownames))))
        r_mat = r.ro.baseenv["colnames<-"](r_mat, r.StrVector(list(map(str, colnames))))
        return r_mat
    
    # # --- Seurat-Specific Converters ---
    # @staticmethod
    # def anndata_to_seurat(adata: AnnData) -> Any:
    #     """
    #     Converts a Python AnnData object to a Seurat object in R memory.
    #     """
    #     r.ro.r('library(Seurat)')
    #     r.ro.r('library(anndata)')
        
    #     r_adata = r.py2r.convert_anndata(adata)
    #     r_seurat_obj = r.ro.r['anndata::as_Seurat'](r_adata)
        
    #     return r_seurat_obj

    # @staticmethod
    # def seurat_to_anndata(r_seurat_obj: Any) -> AnnData:
    #     """
    #     Converts a Seurat object in R back to a Python AnnData object.
    #     """
    #     temp_file = 'temp_seurat_to_anndata.h5ad'
        
    #     r.ro.r('library(anndata)')
        
    #     # We need to assign the R object to a variable in the R environment
    #     r.ro.r.assign("seurat_obj", r_seurat_obj)
    #     r.ro.r(f'anndata::write_h5ad(seurat_obj, filename = "{temp_file}")')
        
    #     adata = ad.read_h5ad(temp_file)
    #     os.remove(temp_file)
        
    #     return adata


# workflow/lib/seurat_converters.py


