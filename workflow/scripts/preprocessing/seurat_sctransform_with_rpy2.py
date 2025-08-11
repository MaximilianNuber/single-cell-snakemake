import os
import sys
import argparse
import anndata as ad
from pyrtools.r_utils import r

# --- Add the 'lib' directory to the Python path ---
script_dir = os.path.dirname(os.path.abspath(__file__))
lib_dir = os.path.join(script_dir, '..', 'lib')
sys.path.insert(0, lib_dir)

# --- Now we can import our utility from the 'lib' directory ---
from ...lib.converters import SeuratRConverters
from ...lib.lazy_r_env import r

# --- Main logic that the script runs ---
def main():
    parser = argparse.ArgumentParser(description="Perform Seurat SCTransform via rpy2.")
    parser.add_argument("--input", required=True, help="Path to the input AnnData H5AD file.")
    parser.add_argument("--output", required=True, help="Path to save the output AnnData object (.h5ad).")
    parser.add_argument("--min-cells", type=int, default=3, help="Minimum number of cells to keep a gene.")
    parser.add_argument("--min-features", type=int, default=200, help="Minimum number of genes to keep a cell.")
    parser.add_argument("--mt-threshold", type=int, default=10, help="Max percent of mitochondrial genes to keep a cell.")
    parser.add_argument("--vst-features", type=int, default=3000, help="Number of features to select with SCTransform.")
    args = parser.parse_args()

    seurat, seurat_object = r.lazy_import_r_packages(["Seurat", "SeuratObject"])

    # Load data in Python
    print(f"Loading AnnData object from: {args.input}")
    adata = ad.read_h5ad(args.input)

    seu = SeuratRConverters.anndata_to_seurat(
        adata,
        min_cells=args.min_cells,
        min_features=args.min_features,
        mt_threshold=args.mt_threshold,
        vst_features=args.vst_features
    )

    col = seurat.PercentageFeatureSet(seu, pattern='^MT-', col.name='percent.mt')
    seu = r.ro.baseenv["[[<-"](seu, "percent.mt", col)

    mask = r.r2py(r.ro.baseenv["[["](seu, "percent.mt"))["percent.mt"].values < args.mt_threshold & r.r2py(r.ro.baseenv["[["](seu, "nCount_RNA"))["nCount_RNA"].values > 0

    r.ro.baseenv["subset"](seu, subset=r.ro.baseenv["percent.mt"] < args.mt_threshold)
    

if __name__ == "__main__":
    main()