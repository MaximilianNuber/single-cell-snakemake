import os
from typing import Any, Optional, Sequence

import pandas as pd
import numpy as np
import anndata 

from .r_env import get_r_environment, r
from .r_converters import RConverters

# A type hint for an AnnData object
AnnData = anndata.AnnData
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