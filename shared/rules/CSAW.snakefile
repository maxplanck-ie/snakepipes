
## CSAW for differential binding / allele-specific binding analysis
rule CSAW:
    input:
        macs2_output = expand("MACS2/{chip_sample}.filtered.BAM_peaks.xls", chip_sample = chip_samples),
        sample_info = sample_info,
        insert_size_metrics =
            "Picard_qc/InsertSizeMetrics/"+chip_samples[0]+".insert_size_metrics.txt" if paired
            else []
    output:
        "CSAW/CSAW.session_info.txt", "CSAW/DiffBinding_analysis.Rdata"
    benchmark:
        "CSAW/.benchmark/CSAW.benchmark"
    params:
        outdir = "CSAW",
        fdr = 0.05,
        paired = paired,
        fragment_length = fragment_length,
        window_size = window_size,
        importfunc = os.path.join(workflow_tools,"snakediff", "R" , "DB_functions.R"),
        allele_info = allele_info
    log: "CSAW/CSAW.log"
    conda: CONDA_ATAC_ENV
    script: "../tools/CSAW.R"
