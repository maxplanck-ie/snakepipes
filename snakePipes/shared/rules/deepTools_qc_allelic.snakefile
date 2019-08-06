### deepTools bamCoverage on allelic BAM files ################################


rule bamCoverage_allelic:
    input:
        bam = "allelic_bams/{sample}.{suffix}.sorted.bam",
        bai = "allelic_bams/{sample}.{suffix}.sorted.bam.bai"
    output:
        "bamCoverage/allele_specific/{sample}.{suffix}.seq_depth_norm.bw"
    params:
        bwBinSize = bwBinSize,
        genome_size = int(genome_size),
        ignoreForNorm = "--ignoreForNormalization " + ignoreForNormalization if ignoreForNormalization else "",
        read_extension = "--extendReads" if paired
                         else "--extendReads " + str(fragmentLength),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
    log:
        out = "bamCoverage/allele_specific/logs/bamCoverage.{sample}.{suffix}.out",
        err = "bamCoverage/allele_specific/logs/bamCoverage.{sample}.{suffix}.err"
    benchmark:
        "bamCoverage/allele_specific/.benchmark/bamCoverage.{sample}.{suffix}.benchmark"
    threads: 16
    conda: CONDA_SHARED_ENV
    shell: bamcov_cmd

### deepTools plotCoverage #####################################################

rule plotCoverage_allelic:
    input:
        bams = expand("allelic_bams/{sample}.{suffix}.sorted.bam", sample=samples, suffix = ['genome1', 'genome2']),
        bais = expand("allelic_bams/{sample}.{suffix}.sorted.bam.bai", sample=samples, suffix = ['genome1', 'genome2'])
    output:
        "deepTools_qc/plotCoverage/read_coverage_allelic.png"
    params:
        labels = " ".join(expand('{sample}.{suffix}', sample=samples, suffix = ['genome1', 'genome2'])),
        read_extension = "--extendReads" if paired
                         else "--extendReads " + str(fragmentLength)
    log:
        out = "deepTools_qc/logs/plotCoverage_allelic.out",
        err = "deepTools_qc/logs/plotCoverage_allelic.err"
    benchmark:
        "deepTools_qc/.benchmark/plotCoverage_allelic.benchmark"
    threads: 24
    conda: CONDA_SHARED_ENV
    shell: plotCoverage_cmd

### deepTools multiBamSummary ##################################################

rule multiBamSummary_allelic:
    input:
        bams = expand("allelic_bams/{sample}.{suffix}.sorted.bam", sample=samples, suffix = ['genome1', 'genome2']),
        bais = expand("allelic_bams/{sample}.{suffix}.sorted.bam.bai", sample=samples, suffix = ['genome1', 'genome2'])
    output:
        "deepTools_qc/multiBamSummary/read_coverage_allelic.bins.npz"
    params:
        labels = " ".join(expand('{sample}.{suffix}', sample=samples, suffix = ['genome1', 'genome2'])),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
        read_extension = "--extendReads" if paired
                         else "--extendReads " + str(fragmentLength)
    log:
        out = "deepTools_qc/logs/multiBamSummary_allelic.out",
        err = "deepTools_qc/logs/multiBamSummary_allelic.err"
    benchmark:
        "deepTools_qc/.benchmark/multiBamSummary_allelic.benchmark"
    threads: 24
    conda: CONDA_SHARED_ENV
    shell: multiBamSummary_cmd

### deepTools plotCorrelation ##################################################

# Pearson: heatmap, scatterplot and correlation matrix
rule plotCorrelation_pearson_allelic:
    input:
        "deepTools_qc/multiBamSummary/read_coverage_allelic.bins.npz"
    output:
        "deepTools_qc/plotCorrelation/correlation.pearson.read_coverage_allelic.tsv"
    params:
        plotcmd = "" if plotFormat == 'None' else
            "--plotFile " + "deepTools_qc/plotCorrelation/correlation.pearson.read_coverage_allelic.heatmap." + plotFormat,
        title='fragment'
    log:
        out = "deepTools_qc/logs/plotCorrelation_pearson_allelic.out",
        err = "deepTools_qc/logs/plotCorrelation_pearson_allelic.err"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_pearson_allelic.benchmark"
    conda: CONDA_SHARED_ENV
    shell: plotCorr_cmd

# Spearman: heatmap, scatterplot and correlation matrix
rule plotCorrelation_spearman_allelic:
    input:
        "deepTools_qc/multiBamSummary/read_coverage_allelic.bins.npz"
    output:
        "deepTools_qc/plotCorrelation/correlation.spearman.read_coverage_allelic.tsv"
    params:
        plotcmd = "" if plotFormat == 'None' else
            "--plotFile " + "deepTools_qc/plotCorrelation/correlation.spearman.read_coverage_allelic.heatmap." + plotFormat,
        title = 'fragment'
    log:
        out = "deepTools_qc/logs/plotCorrelation_spearman_allelic.out",
        err = "deepTools_qc/logs/plotCorrelation_spearman_allelic.err"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_spearman_allelic.benchmark"
    conda: CONDA_SHARED_ENV
    shell: plotCorrSP_cmd

### deepTools plotPCA ##########################################################
rule plotPCA_allelic:
    input:
        "deepTools_qc/multiBamSummary/read_coverage_allelic.bins.npz"
    output:
        "deepTools_qc/plotPCA/PCA.read_coverage_allelic.tsv"
    params:
        plotcmd = "" if plotFormat == 'None' else
                "--plotFile " + "deepTools_qc/plotPCA/PCA.read_coverage_allelic." + plotFormat,
        title='fragment'
    log:
        out = "deepTools_qc/logs/plotPCA_allelic.out",
        err = "deepTools_qc/logs/plotPCA_allelic.err"
    benchmark:
        "deepTools_qc/.benchmark/plotPCA_allelic.benchmark"
    conda: CONDA_SHARED_ENV
    shell: plotPCA_cmd
