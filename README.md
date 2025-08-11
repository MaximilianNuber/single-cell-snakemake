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
