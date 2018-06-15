### deepTools bamCoverage ######################################################


rule bamCoverage:
    input:
        bam = mapping_prg+"/{sample}.bam",
        bai = mapping_prg+"/{sample}.bam.bai"
    output:
        "bamCoverage/{sample}.seq_depth_norm.bw"
    params:
        bw_binsize = bw_binsize,
        genome_size = int(genome_size),
        ignoreForNorm = "--ignoreForNormalization {}".format(ignore_forNorm) if ignore_forNorm else "",
        read_extension = "--extendReads" if paired
                         else "--extendReads {}".format(fragment_length),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
    log:
        out = "bamCoverage/logs/bamCoverage.{sample}.out",
        err = "bamCoverage/logs/bamCoverage.{sample}.err"
    benchmark:
        "bamCoverage/.benchmark/bamCoverage.{sample}.benchmark"
    threads: 16
    conda: CONDA_SHARED_ENV
    shell: bamcov_cmd


### deepTools bamCoverage on filtered BAM files ################################

rule bamCoverage_filtered:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        bai = "filtered_bam/{sample}.filtered.bam.bai"
    output:
        "bamCoverage/{sample}.filtered.seq_depth_norm.bw"
    params:
        bw_binsize = bw_binsize,
        genome_size = int(genome_size),
        ignoreForNorm = "--ignoreForNormalization {}".format(ignore_forNorm) if ignore_forNorm else "",
        read_extension = "--extendReads" if paired
                         else "--extendReads {}".format(fragment_length),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed
                    else "",
    log:
        out = "bamCoverage/logs/bamCoverage.{sample}.filtered.out",
        err = "bamCoverage/logs/bamCoverage.{sample}.filtered.err"
    benchmark:
        "bamCoverage/.benchmark/bamCoverage.{sample}.filtered.benchmark"
    threads: 16
    conda: CONDA_SHARED_ENV
    shell: bamcov_cmd

# TODO: include blacklist!? use deeptools bam filtering options?


### deepTools plotCoverage #####################################################

rule plotCoverage:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample=samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample=samples)
    output:
        "deepTools_qc/plotCoverage/read_coverage.tsv"
    params:
        labels = " ".join(samples),
        read_extension = "--extendReads" if paired
                         else "--extendReads {}".format(fragment_length),
        plotcmd = "" if plot_format == 'None' else
                    "--plotFile deepTools_qc/plotCoverage/read_coverage.{}".format(plot_format)
    log:
        out = "deepTools_qc/logs/plotCoverage.out",
        err = "deepTools_qc/logs/plotCoverage.err"
    benchmark:
        "deepTools_qc/.benchmark/plotCoverage.benchmark"
    threads: 24
    conda: CONDA_SHARED_ENV
    shell: plotCoverage_cmd

### deepTools multiBamSummary ##################################################

rule multiBamSummary:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample=samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample=samples)
    output:
        "deepTools_qc/multiBamSummary/read_coverage.bins.npz"
    params:
        labels = " ".join(samples),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
        read_extension = "--extendReads" if paired
                         else "--extendReads {}".format(fragment_length)
    log:
        out = "deepTools_qc/logs/multiBamSummary.out",
        err = "deepTools_qc/logs/multiBamSummary.err"
    benchmark:
        "deepTools_qc/.benchmark/multiBamSummary.benchmark"
    threads: 24
    conda: CONDA_SHARED_ENV
    shell: multiBamSummary_cmd


### deepTools plotCorrelation ##################################################

# Pearson: heatmap, scatterplot and correlation matrix
rule plotCorrelation_pearson:
    input:
        "deepTools_qc/multiBamSummary/read_coverage.bins.npz"
    output:
        "deepTools_qc/plotCorrelation/correlation.pearson.read_coverage.tsv"
    params:
        plotcmd = "" if plot_format == 'None' else
            "--plotFile deepTools_qc/plotCorrelation/correlation.pearson.read_coverage.heatmap.{}".format(plot_format),
        title='fragment'
    log:
        out = "deepTools_qc/logs/plotCorrelation_pearson.out",
        err = "deepTools_qc/logs/plotCorrelation_pearson.err"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_pearson.benchmark"
    conda: CONDA_SHARED_ENV
    shell: plotCorr_cmd

# Spearman: heatmap, scatterplot and correlation matrix
rule plotCorrelation_spearman:
    input:
        "deepTools_qc/multiBamSummary/read_coverage.bins.npz"
    output:
        "deepTools_qc/plotCorrelation/correlation.spearman.read_coverage.tsv"
    params:
        plotcmd = "" if plot_format == 'None' else
            "--plotFile deepTools_qc/plotCorrelation/correlation.spearman.read_coverage.heatmap.{}".format(plot_format),
        title='fragment'
    log:
        out = "deepTools_qc/logs/plotCorrelation_spearman.out",
        err = "deepTools_qc/logs/plotCorrelation_spearman.err"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_spearman.benchmark"
    conda: CONDA_SHARED_ENV
    shell: plotCorrSP_cmd

### deepTools plotPCA ##########################################################
rule plotPCA:
    input:
        "deepTools_qc/multiBamSummary/read_coverage.bins.npz"
    output:
        "deepTools_qc/plotPCA/PCA.read_coverage.tsv"
    params:
        plotcmd = "" if plot_format == 'None' else
                "--plotFile deepTools_qc/plotPCA/PCA.read_coverage.{}".format(plot_format),
        title='fragment'
    log:
        out = "deepTools_qc/logs/plotPCA.out",
        err = "deepTools_qc/logs/plotPCA.err"
    benchmark:
        "deepTools_qc/.benchmark/plotPCA.benchmark"
    conda: CONDA_SHARED_ENV
    shell: plotPCA_cmd

########## deepTools estimateReadFiltering ###################################

rule estimate_read_filtering:
    input:
        bam = mapping_prg+"/{sample}.bam",
        bai = mapping_prg+"/{sample}.bam.bai"
    output:
        "deepTools_qc/estimateReadFiltering/{sample}_filtering_estimation.txt"
    log:
        out = "deepTools_qc/logs/{sample}.estimateReadFiltering.out",
        err = "deepTools_qc/logs/{sample}.estimateReadFiltering.err"
    conda: CONDA_SHARED_ENV
    shell: estimateReadFiltering_cmd

### deepTools computeGCBias ####################################################

rule computeGCBias:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        bai = "filtered_bam/{sample}.filtered.bam.bai",
    output:
        png = "deepTools_qc/computeGCBias/{sample}.filtered.GCBias.png",
        tsv = "deepTools_qc/computeGCBias/{sample}.filtered.GCBias.freq.tsv"
    params:
        paired = paired,
        fragment_length = fragment_length,
        genome_size = int(genome_size),
        genome_2bit = genome_2bit,
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
        median_fragment_length = "" if paired else "-fragmentLength {}".format(fragment_length),
        sampleSize = downsample if downsample and downsample < 10000000 else 10000000
    log:
        out = "deepTools_qc/logs/computeGCBias.{sample}.filtered.out",
        err = "deepTools_qc/logs/computeGCBias.{sample}.filtered.err"
    benchmark:
        "deepTools_qc/.benchmark/computeGCBias.{sample}.filtered.benchmark"
    threads: 16
    conda: CONDA_SHARED_ENV
    shell: gcbias_cmd

#######InsertSizeMetrics###############
rule bamPE_fragment_size:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample=samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample=samples)
    output:
        "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv"
    params:
        plotcmd = "" if plot_format == 'None' else
                "-o deepTools_qc/bamPEFragmentSize/fragmentSizes.{}".format(plot_format)
    log:
        out = "deepTools_qc/logs/bamPEFragmentSize.out",
        err = "deepTools_qc/logs/bamPEFragmentSize.err"
    threads: 24
    conda: CONDA_SHARED_ENV
    shell: bamPEFragmentSize_cmd
