import os
from typing import Any, Optional, Sequence

import pandas as pd
import numpy as np
import anndata 

from .lazy_r_env import get_r_environment, r

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


import scipy.sparse as sp
class SeuratRConverters:
    """A namespace for Seurat-specific conversion methods."""

    @staticmethod
    def anndata_to_seurat(adata: AnnData, layer: Optional[str] = None, min_cells: int = 0, min_features: int = 0, **kwargs) -> Any:
        """
        Converts a Python AnnData object to a Seurat object in R.
        Handles counts and numeric data correctly.
        """
        seurat, seurat_object = r.lazy_import_r_packages(["Seurat", "SeuratObject"])
        
        # Get the matrix to convert
        mat = adata.X if layer is None else adata.layers[layer]

        # Handle sparse vs. dense matrices
        if sp.issparse(mat):
            # Convert to R sparse matrix
            print("Converting a sparse matrix to Seurat object.")
            r_mat = RConverters.anndata_to_r_matrix(adata.T, layer=layer)
        else:
            # Convert to R dense matrix
            print("Converting a dense matrix to Seurat object.")
            r_mat = RConverters.anndata_to_r_matrix(adata.T, layer=layer)
        
        # Check if matrix contains counts (all integers) or numeric data
        if np.all(np.equal(mat.data if sp.issparse(mat) else mat, (mat.data if sp.issparse(mat) else mat).astype(int))):
            counts = r_mat
            data = r.ro.NULL
            print("Detected counts data.")
        else:
            counts = r.ro.NULL
            data = r_mat
            print("Detected numeric data.")

        assay_obj = seurat_object.CreateAssay5Object(
            counts = counts,
            data = data,
            min_cells = min_cells,
            min_features = min_features,
            csum = r.ro.NULL,
            fsum = r.ro.NULL,
            **kwargs
        )
        obs = r.py2r(adata.obs)
        
        # Create the Seurat object
        r_seurat_obj = seurat_object.CreateSeuratObject(
            counts=assay_obj,
            meta_data = obs
        )
        
        return r_seurat_obj

    @staticmethod
    def seurat_to_anndata(r_seurat_obj: Any) -> AnnData:
        """
        Converts a Seurat object in R back to a Python AnnData object.
        """
        temp_file = 'temp_seurat_to_anndata.h5ad'
        r.ro.r('library(anndata)')
        r.ro.r.assign("seurat_obj", r_seurat_obj)
        r.ro.r(f'anndata::write_h5ad(seurat_obj, filename = "{temp_file}")')
        adata = ad.read_h5ad(temp_file)
        os.remove(temp_file)
        return adata