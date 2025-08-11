# workflow/rules/common.smk

# Define a rule to load the 10x Cell Ranger h5 file and save it
# in a consistent format (e.g., AnnData .h5ad).
# This is a critical starting point for all preprocessing pipelines.
rule load_cellranger_data:
    input:
        lambda wildcards: SAMPLES.loc[wildcards.sample, "cellranger_path"]
    output:
        "results/{sample}/raw_initial_object.h5ad"
    container:
        config['apptainer_image']
    shell:
        """
        python workflow/scripts/common/load_data.py \
            --input {input} \
            --output {output} \
            --sample-id {wildcards.sample}
        """