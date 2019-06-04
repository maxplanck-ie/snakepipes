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

if allele_info == 'FALSE':
    rule calc_matrix_log2r_CSAW_up:
        input:
            csaw_in = "CSAW_{}/CSAW.session_info.txt".format(sample_name),
            bigwigs = expand("deepTools_ChIP/bamCompare/{chip_sample}.filtered.log2ratio.over_{control_name}.bw",zip,chip_sample=chip_samples,control_name=control_names),
            sampleSheet = sampleSheet
        output:
            matrix = "CSAW_{}/CSAW.UP.log2r.matrix".format(sample_name)
        params:
            bigwigs = expand("deepTools_ChIP/bamCompare/{chip_sample}.filtered.log2ratio.over_{control_name}.bw",zip,chip_sample=filtered_dict.keys(),control_name=filtered_dict.values()),
            bed_up = "CSAW_{}/Filtered.results.UP.bed".format(sample_name)
        log:
            out = os.path.join(outdir,"CSAW_{}/logs/deeptools_matrix.log2r.out".format(sample_name)),
            err = os.path.join(outdir,"CSAW_{}/logs/deeptools_matrix.log2r.err".format(sample_name))
        threads: 8
        conda: CONDA_SHARED_ENV
        shell: "if [ -r {params.bed_up} ]; then computeMatrix scale-regions -S {params.bigwigs} -R {params.bed_up} -m 1000 -b 200 -a 200 -o {output.matrix} -p {threads};fi >{log.out} 2>{log.err}"

    rule plot_heatmap_log2r_CSAW_up:
        input:
            matrix = "CSAW_{}/CSAW.UP.log2r.matrix".format(sample_name)
        output:
            image = "CSAW_{}/CSAW.UP.log2r.heatmap.png".format(sample_name),
            sorted_regions = "CSAW_{}/CSAW.UP.log2r.sortedRegions.bed".format(sample_name)
        params:
            smpl_label=' '.join(filtered_dict.keys())
        log:
            out = os.path.join(outdir,"CSAW_{}/logs/deeptools_heatmap.log2r.out".format(sample_name)),
            err = os.path.join(outdir,"CSAW_{}/logs/deeptools_heatmap.log2r.err".format(sample_name))
        conda: CONDA_SHARED_ENV
        shell: "if [ -r {input.matrix} ]; then plotHeatmap --matrixFile {input.matrix} --outFileSortedRegions {output.sorted_regions} --outFileName {output.image} --startLabel Start --endLabel End --legendLocation lower-center -x 'Scaled peak length' --labelRotation 90 --samplesLabel {params.smpl_label} ;fi >{log.out} 2>{log.err}"

    rule calc_matrix_substr_CSAW_up:
        input:
            csaw_in = "CSAW_{}/CSAW.session_info.txt".format(sample_name),
            bigwigs = expand("deepTools_ChIP/bamCompare/{chip_sample}.filtered.subtract.{control_name}.bw",zip,chip_sample=chip_samples,control_name=control_names),
            sampleSheet = sampleSheet
        output:
            matrix = "CSAW_{}/CSAW.UP.substr.matrix".format(sample_name)
        params:
            bigwigs = expand("deepTools_ChIP/bamCompare/{chip_sample}.filtered.subtract.{control_name}.bw",zip,chip_sample=filtered_dict.keys(),control_name=filtered_dict.values()),
            bed_up = "CSAW_{}/Filtered.results.UP.bed".format(sample_name)
        log:
            out = os.path.join(outdir,"CSAW_{}/logs/deeptools_matrix.substr.out".format(sample_name)),
            err = os.path.join(outdir,"CSAW_{}/logs/deeptools_matrix.substr.err".format(sample_name))
        threads: 8
        conda: CONDA_SHARED_ENV
        shell: "if [ -r {params.bed_up} ]; then computeMatrix scale-regions -S {params.bigwigs} -R {params.bed_up} -m 1000 -b 200 -a 200 -o {output.matrix} -p {threads};fi >{log.out} 2>{log.err}"

    rule plot_heatmap_log2r_CSAW_up:
        input:
            matrix = "CSAW_{}/CSAW.UP.substr.matrix".format(sample_name)
        output:
            image = "CSAW_{}/CSAW.UP.substr.heatmap.png".format(sample_name),
            sorted_regions = "CSAW_{}/CSAW.UP.substr.sortedRegions.bed".format(sample_name)
        params:
            smpl_label=' '.join(filtered_dict.keys())
        log:
            out = os.path.join(outdir,"CSAW_{}/logs/deeptools_heatmap.substr.out".format(sample_name)),
            err = os.path.join(outdir,"CSAW_{}/logs/deeptools_heatmap.substr.err".format(sample_name))
        conda: CONDA_SHARED_ENV
        shell: "if [ -r {input.matrix} ]; then plotHeatmap --matrixFile {input.matrix} --outFileSortedRegions {output.sorted_regions} --outFileName {output.image} --startLabel Start --endLabel End --legendLocation lower-center -x 'Scaled peak length' --labelRotation 90 --samplesLabel {params.smpl_label} ;fi >{log.out} 2>{log.err}"
