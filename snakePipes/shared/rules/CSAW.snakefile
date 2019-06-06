sample_name = os.path.splitext(os.path.basename(sampleSheet))[0]
change_direction = ["UP","DOWN","MIXED"]
## CSAW for differential binding / allele-specific binding analysis
rule CSAW:
    input:
        macs2_output = expand("MACS2/{chip_sample}.filtered.BAM_peaks.xls", chip_sample = chip_samples),
        sampleSheet = sampleSheet,
        insert_size_metrics ="deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv" if paired else []
    output:
        "CSAW_{}/CSAW.session_info.txt".format(sample_name), "CSAW_{}/DiffBinding_analysis.Rdata".format(sample_name),expand("CSAW_{}".format(sample_name)+"/Filtered.results.{change_dir}.bed",change_dir=change_direction)
    benchmark:
        "CSAW_{}/.benchmark/CSAW.benchmark".format(sample_name)
    params:
        outdir =os.path.join(outdir,"CSAW_{}".format(sample_name)),
        fdr = fdr,
        absBestLFC=absBestLFC,
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
    rule calc_matrix_log2r_CSAW:
        input:
            csaw_in = "CSAW_{}/CSAW.session_info.txt".format(sample_name),
            bigwigs = expand("deepTools_ChIP/bamCompare/{chip_sample}.filtered.log2ratio.over_{control_name}.bw",zip,chip_sample=filtered_dict.keys(),control_name=filtered_dict.values()),
            sampleSheet = sampleSheet
        output:
            matrix = "CSAW_{}".format(sample_name)+"/CSAW.{change_dir}.log2r.matrix"
        params:
            bed_in = "CSAW_{}".format(sample_name)+"/Filtered.results.{change_dir}.bed"
        log:
            out = os.path.join(outdir,"CSAW_{}".format(sample_name)+"/logs/deeptools_matrix.log2r.{change_dir}.out"),
            err = os.path.join(outdir,"CSAW_{}".format(sample_name)+"/logs/deeptools_matrix.log2r.{change_dir}.err")
        threads: 8
        conda: CONDA_SHARED_ENV
        shell: "if [ -r {params.bed_in} ]; then computeMatrix scale-regions -S {input.bigwigs} -R {params.bed_in} -m 1000 -b 200 -a 200 -o {output.matrix} -p {threads};fi >{log.out} 2>{log.err}"

    rule plot_heatmap_log2r_CSAW:
        input:
            matrix = "CSAW_{}".format(sample_name)+"/CSAW.{change_dir}.log2r.matrix"
        output:
            image = "CSAW_{}".format(sample_name)+"/CSAW.{change_dir}.log2r.heatmap.png",
            sorted_regions = "CSAW_{}".format(sample_name)+"/CSAW.{change_dir}.log2r.sortedRegions.bed"
        params:
            smpl_label=' '.join(filtered_dict.keys())
        log:
            out = os.path.join(outdir,"CSAW_{}".format(sample_name)+"/logs/deeptools_heatmap.log2r.{change_dir}.out"),
            err = os.path.join(outdir,"CSAW_{}".format(sample_name)+"/logs/deeptools_heatmap.log2r{change_dir}.err")
        conda: CONDA_SHARED_ENV
        shell: "if [ -r {input.matrix} ]; then plotHeatmap --matrixFile {input.matrix} --outFileSortedRegions {output.sorted_regions} --outFileName {output.image} --startLabel Start --endLabel End --legendLocation lower-center -x 'Scaled peak length' --labelRotation 90 --samplesLabel {params.smpl_label} ;fi >{log.out} 2>{log.err}"

    rule calc_matrix_cov_CSAW:
        input:
            csaw_in = "CSAW_{}/CSAW.session_info.txt".format(sample_name),
            bigwigs = expand("bamCoverage/{chip_sample}.filtered.seq_depth_norm.bw",chip_sample=filtered_dict.keys()),
            sampleSheet = sampleSheet
        output:
            matrix = "CSAW_{}".format(sample_name)+"/CSAW.{change_dir}.cov.matrix"
        params:
            bed_in = "CSAW_{}".format(sample_name)+"/Filtered.results.{change_dir}.bed"
        log:
            out = os.path.join(outdir,"CSAW_{}".format(sample_name)+"/logs/deeptools_matrix.cov.{change_dir}.out"),
            err = os.path.join(outdir,"CSAW_{}".format(sample_name)+"/logs/deeptools_matrix.cov.{change_dir}.err")
        threads: 8
        conda: CONDA_SHARED_ENV
        shell: "if [ -r {params.bed_in} ]; then computeMatrix scale-regions -S {input.bigwigs} -R {params.bed_in} -m 1000 -b 200 -a 200 -o {output.matrix} -p {threads};fi >{log.out} 2>{log.err}"

    rule plot_heatmap_cov_CSAW:
        input:
            matrix = "CSAW_{}".format(sample_name)+"/CSAW.{change_dir}.cov.matrix"
        output:
            image = "CSAW_{}".format(sample_name)+"/CSAW.{change_dir}.cov.heatmap.png",
            sorted_regions = "CSAW_{}".format(sample_name)+"/CSAW.{change_dir}.cov.sortedRegions.bed"
        params:
            smpl_label=' '.join(filtered_dict.keys())
        log:
            out = os.path.join(outdir,"CSAW_{}".format(sample_name)+"/logs/deeptools_heatmap.cov.{change_dir}.out"),
            err = os.path.join(outdir,"CSAW_{}".format(sample_name)+"/logs/deeptools_heatmap.cov.{change_dir}.err")
        conda: CONDA_SHARED_ENV
        shell: "if [ -r {input.matrix} ]; then plotHeatmap --matrixFile {input.matrix} --outFileSortedRegions {output.sorted_regions} --outFileName {output.image} --startLabel Start --endLabel End --legendLocation lower-center -x 'Scaled peak length' --labelRotation 90 --samplesLabel {params.smpl_label} ;fi >{log.out} 2>{log.err}"
