# workflow/rules/preprocessing.smk

# -----------------------------------------------------------------------------
# Seurat SCTransform Preprocessing Pipeline
# -----------------------------------------------------------------------------
# This rule performs Seurat's SCTransform-based QC and normalization.
# It depends on the output of the 'load_cellranger_data' rule from common.smk.
rule seurat_sctransform:
    input:
        "results/{sample}/raw_initial_object.h5ad"
    output:
        os.path.join(get_workflow_output_dir, "qc_normalized_object.rds")
    params:
        min_cells = config['qc']['min_cells'],
        min_features = config['qc']['min_features'],
        percent_mt_threshold = config['qc']['percent_mt_threshold'],
        vst_n_features = config['seurat_sctransform']['vst_n_features']
    container:
        config['apptainer_image']
    shell:
        """
        Rscript workflow/scripts/preprocessing/seurat_sctransform_script.R \
            --input {input} \
            --output {output} \
            --min-cells {params.min_cells} \
            --min-features {params.min_features} \
            --mt-threshold {params.percent_mt_threshold} \
            --vst-features {params.vst_n_features}
        """

# -----------------------------------------------------------------------------
# Scanpy Preprocessing Pipeline
# -----------------------------------------------------------------------------
# This rule performs Scanpy's standard log-normalization-based QC and normalization.
# It also depends on the output of the 'load_cellranger_data' rule.
rule scanpy_lognorm:
    input:
        "results/{sample}/raw_initial_object.h5ad"
    output:
        os.path.join(get_workflow_output_dir, "qc_normalized_object.h5ad")
    params:
        min_cells = config['qc']['min_cells'],
        min_features = config['qc']['min_features'],
        percent_mt_threshold = config['qc']['percent_mt_threshold'],
        target_sum = config['scanpy_lognorm']['target_sum']
    container:
        config['apptainer_image']
    shell:
        """
        python workflow/scripts/preprocessing/scanpy_lognorm_script.py \
            --input {input} \
            --output {output} \
            --min-cells {params.min_cells} \
            --min-features {params.min_features} \
            --mt-threshold {params.percent_mt_threshold} \
            --target-sum {params.target_sum}
        """