# workflow/rules/scanpy_pipeline.smk

# Include the common rules (e.g., loading data)
include: "common.smk"

# Define the variable that holds the output of the previous step.
# It starts with the output of the common data loading rule.
last_output_file = "results/{sample}/raw_initial_object.h5ad"
# workflow/rules/scanpy_pipeline.smk (updated run_basic_qc rule)

# 1. Basic QC (sc-best-practices style)
rule run_basic_qc:
    input:
        adata = last_output_file,
        raw_matrix = lambda wildcards: SAMPLES.loc[wildcards.sample, "raw_matrix_path"]
    output:
        "results/{sample}/qc_filtered_object.h5ad"
    params:
        min_cells = config['sc_best_practices']['min_cells'],
        min_features = config['sc_best_practices']['min_features'],
    conda:
        "../environments/scanpy.yaml"
    shell:
        """
        python workflow/scripts/pipelines/scanpy_qc.py \
            --adata {input.adata} \
            --matrix {input.raw_matrix} \
            --output {output} \
            --min-cells {params.min_cells} \
            --min-features {params.min_features}
        """

# After the rule, we update the output variable for the next rule.
last_output_file = "results/{sample}/qc_filtered_object.h5ad"
# ... rest of the pipeline ...

# -----------------------------------------------------------------------------
# 2. Optional: EmptyDrops QC
# -----------------------------------------------------------------------------
# This rule is only included if the 'run_empty_drops' flag is True in config.
if config['run_empty_drops']:
    rule run_empty_drops:
        input:
            last_output_file
        output:
            "results/{sample}/empty_drops_filtered_object.h5ad"
        conda:
            "../environments/bioconductor.yaml" # A new environment for Bioconductor
        shell:
            """
            python workflow/scripts/pipelines/run_emptydrops.py \
                --input {input} \
                --output {output}
            """
    # Update the output variable if this rule was included.
    last_output_file = "results/{sample}/empty_drops_filtered_object.h5ad"

# -----------------------------------------------------------------------------
# 3. Normalization
# -----------------------------------------------------------------------------
# This is a required step, but with a choice of methods.
rule run_normalization:
    input:
        last_output_file
    output:
        "results/{sample}/normalized_object.h5ad"
    params:
        method = config['normalization_method']
    conda:
        "../environments/normalization_env.yaml" # An environment with both scanpy and scran
    shell:
        """
        python workflow/scripts/pipelines/run_normalization.py \
            --input {input} \
            --output {output} \
            --method {params.method}
        """

last_output_file = "results/{sample}/normalized_object.h5ad"

# -----------------------------------------------------------------------------
# 4. Final Analysis (PCA, UMAP, Clustering)
# -----------------------------------------------------------------------------
# This is the final step in the pipeline.
rule run_final_scanpy_analysis:
    input:
        last_output_file
    output:
        "results/{sample}/final_clustered_object.h5ad"
    params:
        n_pcs = config['dim_red']['n_pcs'],
        n_neighbors = config['clustering']['n_neighbors'],
        resolution = config['clustering']['resolution']
    conda:
        "../environments/scanpy.yaml"
    shell:
        """
        python workflow/scripts/pipelines/final_scanpy_analysis.py \
            --input {input} \
            --output {output} \
            --n-pcs {params.n_pcs} \
            --n-neighbors {params.n_neighbors} \
            --resolution {params.resolution}
        """

# -----------------------------------------------------------------------------
# 5. The final target rule for this pipeline
# -----------------------------------------------------------------------------
# This rule's input is the output of the final analysis step.
rule all_scanpy:
    input:
        expand("results/{sample}/final_clustered_object.h5ad", sample=SAMPLES.index)