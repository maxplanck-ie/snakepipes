sample_name = os.path.splitext(os.path.basename(sampleSheet))[0]
## CSAW for differential binding / allele-specific binding analysis
rule CSAW:
    input:
        macs2_output = expand("MACS2/{chip_sample}.filtered.BAM_peaks.xls", chip_sample = chip_samples),
        sampleSheet = sampleSheet,
        insert_size_metrics ="deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv" if paired else []
    output:
        "CSAW_{}/CSAW.session_info.txt".format(sample_name), "CSAW_{}/DiffBinding_analysis.Rdata".format(sample_name)
    benchmark:
        "CSAW_{}/.benchmark/CSAW.benchmark".format(sample_name)
    params:
        outdir =os.path.join(outdir,"CSAW_{}".format(sample_name)),
        fdr = 0.05,
        paired = paired,
        fragment_length = fragment_length,
        window_size = window_size,
        importfunc = os.path.join("shared", "rscripts", "DB_functions.R"),
        allele_info = allele_info,
        yaml_path=samples_config,
        insert_size_metrics=os.path.join(outdir,"deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv") if paired else []
    log: 
        out = os.path.join(outdir,"CSAW_{}/logs/CSAW.out".format(sample_name)),
        err = os.path.join(outdir,"CSAW_{}/logs/CSAW.err".format(sample_name))
    conda: CONDA_ATAC_ENV
    script: "../rscripts/CSAW.R"

#rule calc_matrix_CSAW:
#    input:
#        csaw_in = "CSAW_{}/CSAW.session_info.txt".format(sample_name),
#        bams = 
