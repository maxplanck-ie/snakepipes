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

#if not allelic mode
#rule calc_matrix_CSAW_up:
#    input:
#        csaw_in = "CSAW_{}/CSAW.session_info.txt".format(sample_name),
#        bigwigs = expand("deepTools_ChIP/bamCompare/{chip_sample}.filtered.log2ratio.over_{control_name}.bw",chip_sample=chip_samples,control_name=control_names),
#        sampleSheet = sampleSheet
#    output:
#        matrix = "CSAW_{}/CSAW.UP.matrix".format(sample_name)
#        sorted_regions = "CSAW_{}/CSAW.UP.sortedRegions.bed".format(sample_name)
#    params:
#        mdict = dict(zip(chip_sample, control_name)),
#        names_sub = ,
#        bigwigs = lambda wildcards, params:  expand("deepTools_ChIP/bamCompare/{names_sub}.filtered.log2ratio.over_{control_name}.bw",chip_sample=params.names_sub,control_name=params.mdict[params.names_sub]),
#        bed_up = "CSAW_{}/Filtered.results.UP.bed".format(sample_name)
#    log:
#    threads: 8
#    conda: CONDA_SHARED_ENV
#    shell: "if [ -r {params.bed_up}]; then computeMatrix scale-regions -S {params.bigwigs} -R {params.bed_up} -b 1000 -o {output.matrix} --outFileSortedRegions {output.sorted_regions} -p {threads};fi >{log.out} 2>{log.err}"
