### deepTools bamCoverage ######################################################
CONDA_SHARED_ENV = "envs/shared_environment.yaml"

rule bamCoverage:
    input:
        bam = mapping_prg+"/{sample}.bam",
        bai = mapping_prg+"/{sample}.bam.bai"
    output:
        "bamCoverage/{sample}.seq_depth_norm.bw"
    conda:
        CONDA_SHARED_ENV
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
    run:
        shell(bamcov_cmd())


### deepTools bamCoverage on filtered BAM files ################################

rule bamCoverage_filtered:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        bai = "filtered_bam/{sample}.filtered.bam.bai"
    output:
        "bamCoverage/{sample}.filtered.seq_depth_norm.bw"
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
        "bamCoverage/logs/bamCoverage.{sample}.filtered.log"
    benchmark:
        "bamCoverage/.benchmark/bamCoverage.{sample}.filtered.benchmark"
    threads: 16
    run:
        cmd = (bamcov_cmd())
        shell(cmd)

# TODO: include blacklist!? use deeptools bam filtering options?

### deepTools computeGCBias ####################################################

rule computeGCBias:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        bai = "filtered_bam/{sample}.filtered.bam.bai",
        insert_size_metrics =
            "Picard_qc/InsertSizeMetrics/{sample}.insert_size_metrics.txt" if paired
            else []
    output:
        png = "deepTools_qc/computeGCBias/{sample}.filtered.GCBias.png",
        tsv = "deepTools_qc/computeGCBias/{sample}.filtered.GCBias.freq.tsv"
    conda:
        CONDA_SHARED_ENV
    params:
        paired = paired,
        fragment_length = fragment_length,
        genome_size = int(genome_size),
        genome_2bit = genome_2bit,
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else ""
    log:
        "deepTools_qc/logs/computeGCBias.{sample}.filtered.log"
    benchmark:
        "deepTools_qc/.benchmark/computeGCBias.{sample}.filtered.benchmark"
    threads: 16
    run:
        if params.paired:
            median_fragment_length = cf.get_fragment_length(input.insert_size_metrics)
        else:
            median_fragment_length = params.fragment_length

        shell(gcbias_cmd(median_fragment_length))

### deepTools plotCoverage #####################################################

rule plotCoverage:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample=samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample=samples)
    output:
        "deepTools_qc/plotCoverage/read_coverage.png"
    conda:
        CONDA_SHARED_ENV
    params:
        labels = " ".join(samples),
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length)
    log:
        "deepTools_qc/logs/plotCoverage.log"
    benchmark:
        "deepTools_qc/.benchmark/plotCoverage.benchmark"
    threads: 24
    run:
        shell(plotCoverage_cmd())

### deepTools multiBamSummary ##################################################

rule multiBamSummary:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample=samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample=samples)
    output:
        "deepTools_qc/multiBamSummary/read_coverage.bins.npz"
    conda:
        CONDA_SHARED_ENV
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
    run:
        shell(multiBamSummary_cmd())


### deepTools plotCorrelation ##################################################

# Pearson: heatmap, scatterplot and correlation matrix
rule plotCorrelation_pearson:
    input:
        "deepTools_qc/multiBamSummary/read_coverage.bins.npz"
    output:
        heatpng = "deepTools_qc/plotCorrelation/correlation.pearson.read_coverage.heatmap.png",
        tsv = "deepTools_qc/plotCorrelation/correlation.pearson.read_coverage.tsv"
    conda:
        CONDA_SHARED_ENV
    log:
        "deepTools_qc/logs/plotCorrelation_pearson.log"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_pearson.benchmark"
    run:
        shell(plotCorr_cmd('fragment'))

# Spearman: heatmap, scatterplot and correlation matrix
rule plotCorrelation_spearman:
    input:
        "deepTools_qc/multiBamSummary/read_coverage.bins.npz"
    output:
        heatpng = "deepTools_qc/plotCorrelation/correlation.spearman.read_coverage.heatmap.png",
        tsv = "deepTools_qc/plotCorrelation/correlation.spearman.read_coverage.tsv"
    conda:
        CONDA_SHARED_ENV
    log:
        "deepTools_qc/logs/plotCorrelation_spearman.log"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_spearman.benchmark"
    run:
        shell(plotCorrSP_cmd('fragment'))

### deepTools plotPCA ##########################################################
rule plotPCA:
    input:
        "deepTools_qc/multiBamSummary/read_coverage.bins.npz"
    output:
        "deepTools_qc/plotPCA/PCA.read_coverage.png"
    conda:
        CONDA_SHARED_ENV
    log:
        "deepTools_qc/logs/plotPCA.log"
    benchmark:
        "deepTools_qc/.benchmark/plotPCA.benchmark"
    run:
        shell(plotPCA_cmd('fragment'))

########## deepTools estimateReadFiltering ###################################

rule estimate_read_filtering:
    input:
        bam = mapping_prg+"/{sample}.bam",
        bai = mapping_prg+"/{sample}.bam.bai"
    output:
        "deepTools_qc/estimateReadFiltering/{sample}_filtering_estimation.txt"
    conda:
        CONDA_SHARED_ENV
    run:
        shell(estimateReadFiltering_cmd())
