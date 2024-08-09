def get_scaling_factor(sample,input):
    sample_names=[]
    scale_factors=[]
    if os.path.isfile(os.path.join(outdir,input)):
        with open(os.path.join(outdir,input)) as f:
            for idx, line in enumerate(f):
                if idx > 0:
                    sample_names.append(line.split('\t')[0])
                    scale_factors.append((line.split('\t')[1]).rstrip("\n"))
        sf_dict = dict(zip(sample_names, scale_factors))
        scale_factor = sf_dict[sample]

        return float(scale_factor)
    else:
        return float(1)


rule bamCoverage_unique_mappings:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        bai = "filtered_bam/{sample}.filtered.bam.bai"
    output:
        bw_fwd = "bamCoverage/{sample}.uniqueMappings.fwd.bw",
        bw_rev = "bamCoverage/{sample}.uniqueMappings.rev.bw",
    conda:
        CONDA_SHARED_ENV
    params:
        bwBinSize = bwBinSize
    benchmark:
        "bamCoverage/.benchmark/bamCoverage_uniqueMappings.{sample}.benchmark"
    threads: 8
    shell: bamcov_unique_cmd


rule bamCoverage_RPKM:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        bai = "filtered_bam/{sample}.filtered.bam.bai"
    output:
        "bamCoverage/{sample}.RPKM.bw"
    conda:
        CONDA_SHARED_ENV
    params:
        bwBinSize = bwBinSize,
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
        ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else ""
    benchmark:
        "bamCoverage/.benchmark/bamCoverage_RPKM.{sample}.benchmark"
    threads: 8
    shell: bamcov_RPKM_cmd


rule bamCoverage_raw:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        bai = "filtered_bam/{sample}.filtered.bam.bai"
    output:
        "bamCoverage/{sample}.coverage.bw"
    conda:
        CONDA_SHARED_ENV
    params:
        bwBinSize = bwBinSize
    benchmark:
        "bamCoverage/.benchmark/bamCoverage_coverage.{sample}.benchmark"
    threads: 8
    shell: bamcov_raw_cmd


rule multiBamSummary_bed:
    input:
        bam = expand("filtered_bam/{sample}.filtered.bam",sample=samples),
        bai = expand("filtered_bam/{sample}.filtered.bam.bai",sample=samples),
        bed = "Annotation/genes.filtered.bed"
    output:
        scalingFactors = "deepTools_qc/multiBamSummary/scalingFactors.tsv",
        npz = "deepTools_qc/multiBamSummary/results.npz"
    conda:
        CONDA_SHARED_ENV
    params:
        labels = " ".join(samples),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else ""
    benchmark:
        "deepTools_qc/.benchmark/multiBamSummary.bed.benchmark"
    threads: lambda wildcards: 40 if 40<max_thread else max_thread
    shell: multiBamSum_bed_cmd

rule bamCoverage_scaleFactors:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        bai = "filtered_bam/{sample}.filtered.bam.bai",
        scalingFactors = "deepTools_qc/multiBamSummary/scalingFactors.tsv"
    output:
        "bamCoverage/{sample}.scaleFactors.bw"
    conda:
        CONDA_SHARED_ENV
    params:
        scaling_factors = lambda wildcards,input: "--scaleFactor {}".format(get_scaling_factor(wildcards.sample,input.scalingFactors)),
        genome_size = int(genome_size),
        bwBinSize = bwBinSize,
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
        ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
        read_extension = ""
    benchmark:
        "bamCoverage/.benchmark/bamCoverage_scaleFactors.{sample}.benchmark"
    threads: 8
    shell: bamcov_spikein_cmd


rule plotEnrichment:
    input:
        bam = expand("filtered_bam/{sample}.filtered.bam", sample=samples),
        bai = expand("filtered_bam/{sample}.filtered.bam.bai", sample=samples),
        gtf = "Annotation/genes.filtered.gtf"
    output:
        "deepTools_qc/plotEnrichment/plotEnrichment.tsv",
    conda:
        CONDA_SHARED_ENV
    params:
        labels = " ".join(samples),
        plotcmd = "" if plotFormat == 'None' else
            "--plotFile " + "deepTools_qc/plotEnrichment/plotEnrichment." + plotFormat
    benchmark:
        "deepTools_qc/.benchmark/plotEnrichment.benchmark"
    threads: lambda wildcards: 24 if 24<max_thread else max_thread
    shell: plotEnrich_cmd


rule multiBigwigSummary_bed:
    input:
        bw = expand("bamCoverage/{sample}.RPKM.bw", sample=samples),
        bed = "Annotation/genes.filtered.bed",
    output:
        "deepTools_qc/multiBigwigSummary/coverage.bed.npz"
    conda:
        CONDA_SHARED_ENV
    params:
        labels = " ".join(samples)
    benchmark:
        "deepTools_qc/.benchmark/multiBigwigSummary.bed.benchmark"
    threads: 8
    shell: multiBWsum_bed_cmd


# Pearson: heatmap, scatterplot and correlation matrix
rule plotCorr_bed_pearson:
    input:
        "deepTools_qc/multiBigwigSummary/coverage.bed.npz"
    output:
        "deepTools_qc/plotCorrelation/correlation.pearson.bed_coverage.tsv"
    conda:
        CONDA_SHARED_ENV
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_pearson.benchmark"
    params:
        plotcmd = "" if plotFormat == 'None' else
            "--plotFile " + "deepTools_qc/plotCorrelation/correlation.pearson.bed_coverage.heatmap." + plotFormat,
        title='genes'
    shell: plotCorr_cmd


# Spearman: heatmap, scatterplot and correlation matrix
rule plotCorr_bed_spearman:
    input:
        "deepTools_qc/multiBigwigSummary/coverage.bed.npz"
    output:
        "deepTools_qc/plotCorrelation/correlation.spearman.bed_coverage.tsv"
    conda:
        CONDA_SHARED_ENV
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_spearman.benchmark"
    params:
        plotcmd = "" if plotFormat == 'None' else
            "--plotFile " + "deepTools_qc/plotCorrelation/correlation.spearman.bed_coverage.heatmap." + plotFormat,
        title='genes'
    shell: plotCorrSP_cmd


### deepTools plotPCA ##########################################################
rule plotPCA:
    input:
        "deepTools_qc/multiBigwigSummary/coverage.bed.npz"
    output:
        "deepTools_qc/plotPCA/PCA.bed_coverage.tsv"
    conda:
        CONDA_SHARED_ENV
    benchmark:
        "deepTools_qc/.benchmark/plotPCA.benchmark"
    params:
        plotcmd = "" if plotFormat == 'None' else
                "--plotFile " + "deepTools_qc/plotPCA/PCA.bed_coverage." + plotFormat,
        title='genes'
    shell: plotPCA_cmd


########deepTools estimateReadFiltering#########################
rule estimateReadFiltering:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        bai = "filtered_bam/{sample}.filtered.bam.bai",
    output:
        "deepTools_qc/estimateReadFiltering/{sample}_filtering_estimation.txt"
    conda:
        CONDA_SHARED_ENV
    shell: estimateReadFiltering_cmd


#######InsertSizeMetrics###############
rule bamPE_fragment_size:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample=samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample=samples)
    output:
        "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv"
    params:
        plotcmd = "" if plotFormat == 'None' else
                "-o " + "deepTools_qc/bamPEFragmentSize/fragmentSizes." + plotFormat,
    conda:
        CONDA_SHARED_ENV
    threads: lambda wildcards: 24 if 24<max_thread else max_thread
    shell: bamPEFragmentSize_cmd
