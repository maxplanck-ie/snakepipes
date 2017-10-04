
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
    run:
        if params.paired:
            median_fragment_length = cf.get_fragment_length(input.insert_size_metrics)
        else:
            median_fragment_length = params.fragment_length
        shell(
        "( export R_LIBS_USER="+R_libs_path+" && "
        "cat "+os.path.join(workflow_tools,"CSAW.R")+" | "
        ""+os.path.join(R_path,"R")+" --vanilla --slave --args "
        "{input.sample_info} "
        "{params.fdr} "
        "{params.paired} "
        +str(median_fragment_length)+" "
        "{params.window_size} "
        "{params.importfunc} "
        "{params.allele_info} "
        ") 2>&1 | tee {log}")
