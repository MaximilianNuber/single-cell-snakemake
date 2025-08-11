#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(anndataR))
suppressPackageStartupMessages(library(argparser))

# Define and parse arguments
p <- arg_parser("Perform Seurat SCTransform QC and normalization.")
p <- add_argument(p, "--input", help="Path to the input AnnData H5AD file.")
p <- add_argument(p, "--output", help="Path to save the output Seurat object (RDS).")
p <- add_argument(p, "--min-cells", type="integer", default=3, help="Minimum number of cells to keep a gene.")
p <- add_argument(p, "--min-features", type="integer", default=200, help="Minimum number of genes to keep a cell.")
p <- add_argument(p, "--mt-threshold", type="integer", default=10, help="Max percent of mitochondrial genes to keep a cell.")
p <- add_argument(p, "--vst-features", type="integer", default=3000, help="Number of features to select with SCTransform.")
argv <- parse_args(p)

# Load AnnData object from H5AD file
# Note: You'll need the 'anndata' and 'reticulate' R packages installed in your container

create_seurat_object_from_h5ad <- function(
    input_file,
    sample_id = NULL,
    min.cells = 3,
    min.features = 200
) {
    adata <- anndataR::read_h5ad(input_file)
    seurat_obj <- CreateSeuratObject(
        counts = t(adata$X),
        assay = "RNA",
        project = sample_id,
        min.cells = min.cells,
        min.features = min.features
    )
    return(seurat_obj)
}

seurat_obj <- create_seurat_object_from_h5ad(
    input_file = argv$input,
    sample_id = argv$sample_id,
    min.cells = argv$min_cells,
    min.features = argv$min_features
)

# Convert AnnData to a Seurat object
seurat_obj <- CreateSeuratObject(counts = adata$X, project = argv$sample_id, min.cells = argv$min_cells, min.features = argv$min_features)
rm(adata) # Free memory

# Add mitochondrial gene percentage to metadata
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Filter cells based on QC metrics
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > argv$min_features & percent.mt < argv$mt_threshold)

# Run SCTransform
seurat_obj <- SCTransform(seurat_obj, assay = "RNA", verbose = FALSE)

# Save the processed object
saveRDS(seurat_obj, file = argv$output)

print(paste("Seurat object saved to:", argv$output))