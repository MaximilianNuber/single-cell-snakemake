#!/usr/bin/env python3

import argparse
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation
from pathlib import Path

# --- Helper Function for Outlier Detection ---
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

# --- Main Script Logic ---
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--adata", required=True, help="Path to the input h5ad file")
    parser.add_argument("--matrix", required=True, help="Path to the raw matrix file")
    parser.add_argument("--output", required=True, help="Path to the output file")
    parser.add_argument("--min-cells", type=int, default=3, help="Min cells for gene filtering.")
    parser.add_argument("--min-features", type=int, default=200, help="Min features for cell filtering.")
    args = parser.parse_args()

    input_path = Path(args.adata)
    matrix_path = Path(args.matrix)
    output_path = Path(args.output)
    
    # Derive sample name from input filename (e.g. data/S01.h5ad -> S01)
    sample = input_path.stem
    
    print("Running QC for:", sample)
    print("  Adata:  ", input_path)
    print("  Matrix: ", matrix_path)
    print("  Output: ", output_path)

    try:
        adata = sc.read_h5ad(input_path)
        
        # Add gene information
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
        adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
        
        # Calculate QC metrics
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
        )
        
        # Filter based on statistical outliers (using MAD)
        adata.obs["outlier"] = (
            is_outlier(adata, "log1p_total_counts", 5)
            | is_outlier(adata, "log1p_n_genes_by_counts", 5)
            | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
        )
        
        adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
            adata.obs["pct_counts_mt"] > 8
        )
        
        adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()

        # Final filtering based on min genes/cells
        # sc.pp.filter_cells(adata, min_genes=args.min_features)
        sc.pp.filter_genes(adata, min_cells=args.min_cells)

        # Save the filtered object
        print(f"Saving QC-filtered object to: {output_path}")
        adata.write(output_path, compression="gzip")
        
    except Exception as e:
        print(f"An error occurred during QC for sample {sample}: {e}")
        # Re-raise the exception to fail the Snakemake job
        raise

if __name__ == "__main__":
    main()