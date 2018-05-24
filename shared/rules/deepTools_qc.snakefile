### deepTools bamCoverage ######################################################
CONDA_SHARED_ENV = "envs/shared_environment.yaml"

rule bamCoverage:
    input:
        bam = mapping_prg+"/{sample}.bam",
        bai = mapping_prg+"/{sample}.bam.bai"
    output:
        "bamCoverage/{sample}.seq_depth_norm.bw"
    params:
        bw_binsize = bw_binsize,
        genome_size = int(genome_size),
        ignoreForNorm = "--ignoreForNormalization " + ignore_forNorm if ignore_forNorm else "",
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length)
    log:
        "bamCoverage/logs/bamCoverage.{sample}.log"
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
        ignoreForNorm = "--ignoreForNormalization " + ignore_forNorm if ignore_forNorm else "",
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
    log:
        "bamCoverage/logs/bamCoverage.{sample}.filtered.log"
    benchmark:
        "bamCoverage/.benchmark/bamCoverage.{sample}.filtered.benchmark"
    threads: 16
    conda: CONDA_SHARED_ENV
    shell: bamcov_cmd

# TODO: include blacklist!? use deeptools bam filtering options?

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
        blacklist = "--blackListFileName "+ blacklist_bed if blacklist_bed
                    else "",
        median_fragment_length = "" if paired else "-fragmentLength " + fragment_length
    log:
        "deepTools_qc/logs/computeGCBias.{sample}.filtered.log"
    benchmark:
        "deepTools_qc/.benchmark/computeGCBias.{sample}.filtered.benchmark"
    threads: 16
    conda: CONDA_SHARED_ENV
    shell: gcbias_cmd


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
                         else "--extendReads "+str(fragment_length),
        plotcmd = "" if plot_format == 'None' else
                    "--plotFile " + "deepTools_qc/plotCoverage/read_coverage." + plot_format
    log:
        "deepTools_qc/logs/plotCoverage.log"
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
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length)
    log:
        "deepTools_qc/logs/multiBamSummary.log"
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
            "--plotFile " + "deepTools_qc/plotCorrelation/correlation.pearson.read_coverage.heatmap." + plot_format,
        title='fragment'
    log:
        "deepTools_qc/logs/plotCorrelation_pearson.log"
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
            "--plotFile " + "deepTools_qc/plotCorrelation/correlation.spearman.read_coverage.heatmap." + plot_format,
        title='fragment'
    log:
        "deepTools_qc/logs/plotCorrelation_spearman.log"
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
                "--plotFile " + "deepTools_qc/plotPCA/PCA.read_coverage." + plot_format,
        title='fragment'
    log:
        "deepTools_qc/logs/plotPCA.log"
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
    conda: CONDA_SHARED_ENV
    shell: estimateReadFiltering_cmd

#######InsertSizeMetrics###############
rule bamPE_fragment_size:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample=samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample=samples)
    output:
        "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv"
    params:
        plotcmd = "" if plot_format == 'None' else
                "--plotFile " + "deepTools_qc/bamPEFragmentSize/fragmentSizes." + plot_format,
    log:
        "deepTools_qc/logs/bamPEFragmentSize.log"
    threads: 24
    conda: CONDA_SHARED_ENV
    shell: bamPEFragmentSize_cmd
