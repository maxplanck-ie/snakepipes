

rule bamCoverage_RPKM:
    input:
        bam = mapping_prg+"/{sample}.bam",
        bai = mapping_prg+"/{sample}.bam.bai"
    output:
        "bamCoverage/{sample}.RPKM.bw"
    conda:
        CONDA_SHARED_ENV
    params:
        bw_binsize = bw_binsize
    log:
        "bamCoverage/logs/bamCoverage_RPKM.{sample}.log"
    benchmark:
        "bamCoverage/.benchmark/bamCoverage_RPKM.{sample}.benchmark"
    threads: 8
    shell: bamcov_RPKM_cmd


rule bamCoverage_raw:
    input:
        bam = mapping_prg+"/{sample}.bam",
        bai = mapping_prg+"/{sample}.bam.bai"
    output:
        "bamCoverage/{sample}.coverage.bw"
    conda:
        CONDA_SHARED_ENV
    params:
        bw_binsize = bw_binsize
    log:
        "bamCoverage/logs/bamCoverage_coverage.{sample}.log"
    benchmark:
        "bamCoverage/.benchmark/bamCoverage_coverage.{sample}.benchmark"
    threads: 8
    shell: bamcov_raw_cmd


rule plotEnrichment:
    input:
        bam = expand(mapping_prg+"/{sample}.bam", sample=samples),
        bai = expand(mapping_prg+"/{sample}.bam.bai", sample=samples),
        gtf = "Annotation/genes.filtered.gtf",
        gtf2= "Annotation/genes.filtered.transcripts.gtf"
    output:
        png = "deepTools_qc/plotEnrichment/plotEnrichment.png",
        tsv = "deepTools_qc/plotEnrichment/plotEnrichment.tsv",
    conda:
        CONDA_SHARED_ENV
    params:
        labels = " ".join(samples),
    log:
        "deepTools_qc/logs/plotEnrichment.log"
    benchmark:
        "deepTools_qc/.benchmark/plotEnrichment.benchmark"
    threads: 8
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
    log:
        "deepTools_qc/logs/multiBigwigSummary.log"
    benchmark:
        "deepTools_qc/.benchmark/multiBigwigSummary.bed.benchmark"
    threads: 8
    shell: multiBWsum_bed_cmd


# Pearson: heatmap, scatterplot and correlation matrix
rule plotCorr_bed_pearson:
    input:
        "deepTools_qc/multiBigwigSummary/coverage.bed.npz"
    output:
        heatpng = "deepTools_qc/plotCorrelation/correlation.pearson.bed_coverage.heatmap.png",
        #scatterpng = "deepTools_qc/plotCorrelation/correlation.pearson.bed_coverage.scatterplot.png",
        tsv = "deepTools_qc/plotCorrelation/correlation.pearson.bed_coverage.tsv"
    conda:
        CONDA_SHARED_ENV
    log:
        "deepTools_qc/logs/plotCorrelation_pearson.log"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_pearson.benchmark"
    params: label='gene'
    shell: plotCorr_cmd


# Spearman: heatmap, scatterplot and correlation matrix
rule plotCorr_bed_spearman:
    input:
        "deepTools_qc/multiBigwigSummary/coverage.bed.npz"
    output:
        heatpng = "deepTools_qc/plotCorrelation/correlation.spearman.bed_coverage.heatmap.png",
        #scatterpng = "deepTools_qc/plotCorrelation/correlation.spearman.bed_coverage.scatterplot.png",
        tsv = "deepTools_qc/plotCorrelation/correlation.spearman.bed_coverage.tsv"
    conda:
        CONDA_SHARED_ENV
    log:
        "deepTools_qc/logs/plotCorrelation_spearman.log"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_spearman.benchmark"
    params: label='gene'
    shell: plotCorrSP_cmd


### deepTools plotPCA ##########################################################
rule plotPCA:
    input:
        "deepTools_qc/multiBigwigSummary/coverage.bed.npz"
    output:
        "deepTools_qc/plotPCA/PCA.bed_coverage.png"
    conda:
        CONDA_SHARED_ENV
    log:
        "deepTools_qc/logs/plotPCA.log"
    benchmark:
        "deepTools_qc/.benchmark/plotPCA.benchmark"
    params: label='gene'
    shell: plotPCA_cmd


########deepTools estimateReadFiltering#########################
rule estimateReadFiltering:
    input:
        bam = mapping_prg+"/{sample}.bam",
        bai = mapping_prg+"/{sample}.bam.bai",
    output:
        "deepTools_qc/estimateReadFiltering/{sample}_filtering_estimation.txt"
    conda:
        CONDA_SHARED_ENV
    shell: estimateReadFiltering_cmd


#######InsertSizeMetrics###############
rule bamPE_fragment_size:
    input:
        bams = expand(mapping_prg+"/{sample}.bam", sample=samples),
        bais = expand(mapping_prg+"/{sample}.bam.bai", sample=samples)
    output:
        "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv"
    conda:
        CONDA_SHARED_ENV
    log:
        "deepTools_qc/bamPEFragmentSize/log"
    threads: 24
    shell: bamPEFragmentSize_cmd
