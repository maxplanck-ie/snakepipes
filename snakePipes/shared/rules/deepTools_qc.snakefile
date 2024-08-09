### deepTools bamCoverage ######################################################


rule bamCoverage:
    input:
        bam = aligner+"/{sample}.bam",
        bai = aligner+"/{sample}.bam.bai"
    output:
        "bamCoverage/{sample}.seq_depth_norm.bw"
    params:
        bwBinSize = bwBinSize,
        genome_size = int(genome_size),
        ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads {}".format(fragmentLength),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
        scaling_factors = ""
    benchmark:
        "bamCoverage/.benchmark/bamCoverage.{sample}.benchmark"
    threads: lambda wildcards: 16 if 16<max_thread else max_thread  # 4GB per core
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
        bwBinSize = bwBinSize,
        genome_size = int(genome_size),
        ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads {}".format(fragmentLength),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed
                    else "",
        scaling_factors ="",
        binSize = ""
    benchmark:
        "bamCoverage/.benchmark/bamCoverage.{sample}.filtered.benchmark"
    threads: lambda wildcards: 16 if 16<max_thread else max_thread  # 4GB per core
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
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads {}".format(fragmentLength),
        plotcmd = "" if plotFormat == 'None' else
                    "--plotFile deepTools_qc/plotCoverage/read_coverage.{}".format(plotFormat)
    benchmark:
        "deepTools_qc/.benchmark/plotCoverage.benchmark"
    threads: lambda wildcards: 24 if 24<max_thread else max_thread
    conda: CONDA_SHARED_ENV
    shell: plotCoverage_cmd

### deepTools multiBamSummary ##################################################

rule multiBamSummary:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample=samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample=samples)
    output:
        npz = "deepTools_qc/multiBamSummary/read_coverage.bins.npz"
    params:
        labels = " ".join(samples),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads {}".format(fragmentLength),
        scaling_factors = "--scalingFactors deepTools_qc/multiBamSummary/scaling_factors.txt",
        binSize = "",
        spikein_region = ""
    benchmark:
        "deepTools_qc/.benchmark/multiBamSummary.benchmark"
    threads: lambda wildcards: 40 if 40<max_thread else max_thread
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
        plotcmd = "" if plotFormat == 'None' else
            "--plotFile deepTools_qc/plotCorrelation/correlation.pearson.read_coverage.heatmap.{}".format(plotFormat),
        title='fragment'
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
        plotcmd = "" if plotFormat == 'None' else
            "--plotFile deepTools_qc/plotCorrelation/correlation.spearman.read_coverage.heatmap.{}".format(plotFormat),
        title='fragment'
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
        plotcmd = "" if plotFormat == 'None' else
                "--plotFile deepTools_qc/plotPCA/PCA.read_coverage.{}".format(plotFormat),
        title='fragment'
    benchmark:
        "deepTools_qc/.benchmark/plotPCA.benchmark"
    conda: CONDA_SHARED_ENV
    shell: plotPCA_cmd

########## deepTools estimateReadFiltering ###################################

rule estimate_read_filtering:
    input:
        bam = aligner+"/{sample}.bam",
        bai = aligner+"/{sample}.bam.bai"
    output:
        "deepTools_qc/estimateReadFiltering/{sample}_filtering_estimation.txt"
    conda: CONDA_SHARED_ENV
    shell: estimateReadFiltering_cmd

### deepTools computeGCBias ####################################################

rule computeGCBias:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        bai = "filtered_bam/{sample}.filtered.bam.bai"
    output:
        png = "deepTools_qc/computeGCBias/{sample}.filtered.GCBias.png",
        tsv = "deepTools_qc/computeGCBias/{sample}.filtered.GCBias.freq.tsv"
    params:
        pairedEnd = pairedEnd,
        fragmentLength = fragmentLength,
        genome_size = int(genome_size),
        genome_2bit = genome_2bit,
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
        median_fragment_length = "" if pairedEnd else "--fragmentLength {}".format(fragmentLength),
        sampleSize = downsample if downsample and downsample < 10000000 else 10000000
    benchmark:
        "deepTools_qc/.benchmark/computeGCBias.{sample}.filtered.benchmark"
    threads: lambda wildcards: 16 if 16<max_thread else max_thread
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
        plotcmd = "" if plotFormat == 'None' else
                "-o deepTools_qc/bamPEFragmentSize/fragmentSizes.{}".format(plotFormat)
    threads: lambda wildcards: 24 if 24<max_thread else max_thread
    conda: CONDA_SHARED_ENV
    shell: bamPEFragmentSize_cmd
