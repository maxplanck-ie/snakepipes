### deepTools bamCoverage on allelic BAM files ################################
CONDA_SHARED_ENV = "envs/shared_environment.yaml"

rule bamCoverage_allelic:
    input:
        bam = "allelic_bams/{sample}.{suffix}.sorted.bam",
        bai = "allelic_bams/{sample}.{suffix}.sorted.bam.bai"
    output:
        "bamCoverage/allele_specific/{sample}.{suffix}.seq_depth_norm.bw"
    conda:
        CONDA_SHARED_ENV
    params:
        bw_binsize = bw_binsize,
        genome_size = int(genome_size),
        ignoreForNorm = "--ignoreForNormalization " + ignore_forNorm if ignore_forNorm else "",
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
    log:
        "bamCoverage/allele_specific/logs/bamCoverage.{sample}.{suffix}.log"
    benchmark:
        "bamCoverage/allele_specific/.benchmark/bamCoverage.{sample}.{suffix}.benchmark"
    threads: 16
    shell: bamcov_cmd + " {params.blacklist}"


### deepTools computeGCBias ####################################################

rule computeGCBias_allelic:
    input:
        bam = "allelic_bams/{sample}.{suffix}.sorted.bam",
        bai = "allelic_bams/{sample}.{suffix}.sorted.bam.bai",
        insert_size_metrics =
            "Picard_qc/InsertSizeMetrics/{sample}.insert_size_metrics.txt" if paired
            else []
    output:
        png = "deepTools_qc/computeGCBias/{sample}.{suffix}.GCBias.png",
        tsv = "deepTools_qc/computeGCBias/{sample}.{suffix}.GCBias.freq.tsv"
    conda:
        CONDA_SHARED_ENV
    params:
        paired = paired,
        fragment_length = fragment_length,
        genome_size = int(genome_size),
        genome_2bit = genome_2bit,
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else ""
        median_fragment_length = cf.get_fragment_length(input.insert_size_metrics) if paired else fragment_length
    log:
        "deepTools_qc/logs/computeGCBias.{sample}.{suffix}.log"
    benchmark:
        "deepTools_qc/.benchmark/computeGCBias.{sample}.{suffix}.benchmark"
    threads: 16
    params: median_fragment_length=cf.get_fragment_length(input.insert_size_metrics) if paired else params.fragment_length
    shell: gcbias_cmd

### deepTools plotCoverage #####################################################

rule plotCoverage_allelic:
    input:
        bams = expand("allelic_bams/{sample}.{suffix}.sorted.bam", sample=samples, suffix = ['genome1', 'genome2']),
        bais = expand("allelic_bams/{sample}.{suffix}.sorted.bam.bai", sample=samples, suffix = ['genome1', 'genome2'])
    output:
        "deepTools_qc/plotCoverage/read_coverage_allelic.png"
    conda:
        CONDA_SHARED_ENV
    params:
        labels = " ".join(samples),
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length)
    log:
        "deepTools_qc/logs/plotCoverage_allelic.log"
    benchmark:
        "deepTools_qc/.benchmark/plotCoverage_allelic.benchmark"
    threads: 24
    shell: plotCoverage_cmd

### deepTools multiBamSummary ##################################################

rule multiBamSummary_allelic:
    input:
        bams = expand("allelic_bams/{sample}.{suffix}.sorted.bam", sample=samples, suffix = ['genome1', 'genome2']),
        bais = expand("allelic_bams/{sample}.{suffix}.sorted.bam.bai", sample=samples, suffix = ['genome1', 'genome2'])
    output:
        "deepTools_qc/multiBamSummary/read_coverage_allelic.bins.npz"
    conda:
        CONDA_SHARED_ENV
    params:
        labels = " ".join(expand('{sample}.{suffix}', sample=samples, suffix = ['genome1', 'genome2']) ),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length)
    log:
        "deepTools_qc/logs/multiBamSummary_allelic.log"
    benchmark:
        "deepTools_qc/.benchmark/multiBamSummary_allelic.benchmark"
    threads: 24
    shell: multiBamSummary_cmd

### deepTools plotCorrelation ##################################################

# Pearson: heatmap, scatterplot and correlation matrix
rule plotCorrelation_pearson_allelic:
    input:
        "deepTools_qc/multiBamSummary/read_coverage_allelic.bins.npz"
    output:
        heatpng = "deepTools_qc/plotCorrelation/correlation.pearson.read_coverage_allelic.heatmap.png",
        tsv = "deepTools_qc/plotCorrelation/correlation.pearson.read_coverage_allelic.tsv"
    conda:
        CONDA_SHARED_ENV
    log:
        "deepTools_qc/logs/plotCorrelation_pearson_allelic.log"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_pearson_allelic.benchmark"
    params: label='fragment'
    shell: plotCorr_cmd

# Spearman: heatmap, scatterplot and correlation matrix
rule plotCorrelation_spearman_allelic:
    input:
        "deepTools_qc/multiBamSummary/read_coverage_allelic.bins.npz"
    output:
        heatpng = "deepTools_qc/plotCorrelation/correlation.spearman.read_coverage_allelic.heatmap.png",
        tsv = "deepTools_qc/plotCorrelation/correlation.spearman.read_coverage_allelic.tsv"
    conda:
        CONDA_SHARED_ENV
    log:
        "deepTools_qc/logs/plotCorrelation_spearman_allelic.log"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_spearman_allelic.benchmark"
    params: label='fragment'
    shell: plotCorrSP_cmd

### deepTools plotPCA ##########################################################
rule plotPCA_allelic:
    input:
        "deepTools_qc/multiBamSummary/read_coverage_allelic.bins.npz"
    output:
        "deepTools_qc/plotPCA/PCA.read_coverage_allelic.png"
    conda:
        CONDA_SHARED_ENV
    log:
        "deepTools_qc/logs/plotPCA_allelic.log"
    benchmark:
        "deepTools_qc/.benchmark/plotPCA_allelic.benchmark"
    params: label='fragment'
    shell: plotPCA_cmd
