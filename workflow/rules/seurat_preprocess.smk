rule runSamplePreprocessing:
    '''
    This rule preprcess the individual samples from the cellranger output.
    '''
    input:
        # result_folder = lambda wildcards: SAMPLES[(wildcards.sample_name)]["out_cellranger"]
        result_folder = rules.runCellrangerMultiRun.output.folder
    output:
        rds = config["out_location"] + "Seurat/object/{sample_name}_obj_preQC.rds",
        meta = config["out_location"] + "Seurat/table/{sample_name}_meta_preQC.tsv",
    conda: config["env_seurat"]
    log:
        'logs/Seurat/{sample_name}/runSamplePreprocessing.log'
    benchmark:
        'benchmarks/Seurat/{sample_name}/runSamplePreprocessing.txt'
    params:
        id_org = config["ref"]["organism_id"]
    script:
        "../scripts/01_sample_preprocessing_snakemake.R"

rule runAggregateQC:
    '''
    This rule aggreagates the metadata form the individual samples and generates the plot to define the QC thresholds.
    '''
    input:
        # this is needed to trigger it after the generation of the outputs on all the outputs
        meta_tables = expand(rules.runSamplePreprocessing.output.meta,sample_name=SAMPLES_merge.keys())
    output:
        plot_mito = config["out_location"] + "Seurat/plot/fixed_histo_mito_V5.pdf",
        plot_feature = config["out_location"] + "Seurat/plot/fixed_histo_features_V5.pdf",
        LUT_QC = config["out_location"] + "Seurat/table/LUT_QC.csv",
        meta_total = config["out_location"] + "Seurat/table/meta_total_beforeQC_V5.tsv"
    conda: config["env_seurat"]
    log:
        'logs/Seurat/runAggregateQC.log'
    benchmark:
        'benchmarks/Seurat/runAggregateQC.txt'
    params:
        # id_org = config["ref"]["organism_id"]
    script:
        "../scripts/01_aggregate_QC_snakemake.R"

rule runApplyQC:
    '''
    This rule preprcess the individual samples from the cellranger output.
    '''
    input:
        rds = rules.runSamplePreprocessing.output.rds,
        LUT_QC = rules.runAggregateQC.output.LUT_QC
    output:
        rds = config["out_location"] + "Seurat/object/{sample_name}_obj_postQC.rds",
        meta = config["out_location"] + "Seurat/table/{sample_name}_meta_postQC.tsv",
    conda: config["env_seurat"]
    log:
        'logs/Seurat/{sample_name}/runApplyQC.log'
    benchmark:
        'benchmarks/Seurat/{sample_name}/runApplyQC.txt'
    params:
        # id_org = config["ref"]["organism_id"]
    script:
        "../scripts/02_apply_QC_snakemake.R"
