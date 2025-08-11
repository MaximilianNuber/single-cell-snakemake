# Define the project root directory name
PROJECT_ROOT="snakemake-sc-pipeline"

echo "Creating project directory structure in: $(pwd)/$PROJECT_ROOT"

# Create main project directory
mkdir -p "$PROJECT_ROOT"
cd "$PROJECT_ROOT"

# Create core directories
mkdir -p config workflow/rules workflow/scripts/common workflow/scripts/preprocessing workflow/scripts/downstream workflow/environments apptainer results

# Create .gitignore
cat << 'EOF' > .gitignore
# Ignore results and large data files
results/
.snakemake/
.apptainer/
*.sif
*.log
*.tmp
*.temp
*.bak
*.orig
EOF

# Create LICENSE (MIT License as a common example)
cat << 'EOF' > LICENSE
MIT License

Copyright (c) 2025 Your Name or Organization

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
EOF

# Create README.md (initial placeholder)
cat << 'EOF' > README.md
# Single-Cell Preprocessing Snakemake Workflow

This repository contains a general Snakemake workflow for single-cell RNA sequencing preprocessing, designed with "Workflows as Software" principles in mind. It offers flexible options for different preprocessing pipelines (e.g., Seurat SCTransform, Scanpy).

## Getting Started

1.  **Clone this repository:**
    ```bash
    git clone [https://github.com/your_username/snakemake-sc-pipeline.git](https://github.com/your_username/snakemake-sc-pipeline.git)
    cd snakemake-sc-pipeline
    ```

2.  **Configure your workflow:**
    * Edit `config/config.yaml` to specify your desired preprocessing method and parameters.
    * Edit `config/samples.tsv` to list your input 10x Cell Ranger output paths.

3.  **Build the Apptainer container:**
    Follow the instructions in the `apptainer/` directory to build the necessary container image.

4.  **Run the workflow:**
    (Instructions will be added here later)

## Workflow Structure

* `config/`: Configuration files (`config.yaml`, `samples.tsv`)
* `workflow/`: Snakemake rules, scripts, and environments
* `apptainer/`: Apptainer definition files
* `results/`: Output directory (ignored by Git)

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.
EOF

# Create config/config.yaml
cat << 'EOF' > config/config.yaml
# --- Workflow Options ---
# Select the preprocessing pipeline to run.
# Available options: "seurat_sctransform", "scanpy_lognorm"
preprocessing_workflow: "seurat_sctransform"

# --- QC and Filtering Parameters (shared across workflows) ---
qc:
  min_cells: 3               # Filter genes expressed in fewer than N cells
  min_features: 200          # Filter cells with fewer than N genes
  percent_mt_threshold: 10   # Filter cells with more than N% mitochondrial genes

# --- Preprocessing-specific Parameters ---
seurat_sctransform:
  vst_n_features: 3000       # Number of features to select for SCTransform

scanpy_lognorm:
  target_sum: 1e4            # Total count to normalize each cell to

# --- Downstream Analysis Parameters (shared) ---
clustering:
  resolution: 0.8
  n_neighbors: 30
  n_pcs: 50
  
# --- Software and Container Paths ---
# Path to the pre-built Apptainer container
apptainer_image: "apptainer/sc_pipeline.sif"
EOF

# Create config/samples.tsv
cat << 'EOF' > config/samples.tsv
sample_id	cellranger_path
sample_1	/path/to/my_data/cellranger_output_1/filtered_feature_bc_matrix.h5
sample_2	/path/to/my_data/cellranger_output_2/filtered_feature_bc_matrix.h5
EOF

echo "Directory structure and initial config files created successfully in $PROJECT_ROOT."
echo "You are now inside the '$PROJECT_ROOT' directory."
echo "Next, we'll work on the workflow/Snakefile."
