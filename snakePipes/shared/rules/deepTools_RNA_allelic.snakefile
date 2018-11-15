

rule bamCoverage_RPKM_allelic:
    input:
        bam = "allelic_bams/{sample}.{suffix}.sorted.bam",
        bai = "allelic_bams/{sample}.{suffix}.sorted.bam.bai"
    output:
        "bamCoverage/allele_specific/{sample}.{suffix}.RPKM.bw"
    conda:
        CONDA_SHARED_ENV
    params:
        bw_binsize = config["bw_binsize"]
    log:
        out="bamCoverage/allele_specific/logs/bamCoverage_RPKM.{sample}.{suffix}.out",
        err="bamCoverage/allele_specific/logs/bamCoverage_RPKM.{sample}.{suffix}.err",
    benchmark:
        "bamCoverage/allele_specific/.benchmark/bamCoverage_RPKM.{sample}.{suffix}.benchmark"
    threads: 8
    shell: bamcov_RPKM_cmd


rule bamCoverage_raw_allelic:
    input:
        bam = "allelic_bams/{sample}.{suffix}.sorted.bam",
        bai = "allelic_bams/{sample}.{suffix}.sorted.bam.bai"
    output:
        "bamCoverage/allele_specific/{sample}.{suffix}.coverage.bw"
    conda:
        CONDA_SHARED_ENV
    params:
        bw_binsize = bw_binsize
    log:
        out="bamCoverage/allele_specificlogs/bamCoverage_coverage.{sample}.{suffix}.out",
        err="bamCoverage/allele_specificlogs/bamCoverage_coverage.{sample}.{suffix}.err"
    benchmark:
        "bamCoverage/allele_specific.benchmark/bamCoverage_coverage.{sample}.{suffix}.benchmark"
    threads: 8
    shell: bamcov_raw_cmd


rule plotEnrichment_allelic:
    input:
        bam = expand("allelic_bams/{sample}.{suffix}.sorted.bam", sample=samples, suffix = ['genome1', 'genome2']),
        bai = expand("allelic_bams/{sample}.{suffix}.sorted.bam.bai", sample=samples, suffix = ['genome1', 'genome2']),
        gtf = "Annotation/genes.filtered.gtf"
    output:
        "deepTools_qc/plotEnrichment/plotEnrichment_allelic.tsv"
    conda:
        CONDA_SHARED_ENV
    params:
        labels = " ".join(expand('{sample}.{suffix}', sample=samples, suffix = ['genome1', 'genome2'])),
        plotcmd = "" if plot_format == 'None' else
            "--plotFile " + "deepTools_qc/plotEnrichment/plotEnrichment_allelic." + plot_format
    log:
        out="deepTools_qc/logs/plotEnrichment_allelic.out",
        err="deepTools_qc/logs/plotEnrichment_allelic.err"
    benchmark:
        "deepTools_qc/.benchmark/plotEnrichment_allelic.benchmark"
    threads: 24
    shell: plotEnrich_cmd


rule multiBigwigSummary_bed_allelic:
    input:
        bw = expand("bamCoverage/allele_specific/{sample}.{suffix}.RPKM.bw", sample=samples, suffix = ['genome1', 'genome2']),
        bed = "Annotation/genes.filtered.bed"
    output:
        "deepTools_qc/multiBigwigSummary/coverage_allelic.bed.npz"
    conda:
        CONDA_SHARED_ENV
    params:
        labels = " ".join(expand('{sample}.{suffix}', sample=samples, suffix = ['genome1', 'genome2']))
    log:
        out="deepTools_qc/logs/multiBigwigSummary_allelic.out",
        err="deepTools_qc/logs/multiBigwigSummary_allelic.err"
    benchmark:
        "deepTools_qc/.benchmark/multiBigwigSummary.bed_allelic.benchmark"
    threads: 8
    shell: multiBWsum_bed_cmd


# Pearson: heatmap, scatterplot and correlation matrix
rule plotCorr_bed_pearson_allelic:
    input:
        "deepTools_qc/multiBigwigSummary/coverage_allelic.bed.npz"
    output:
        "deepTools_qc/plotCorrelation/correlation.pearson.bed_coverage_allelic.tsv"
    conda:
        CONDA_SHARED_ENV
    log:
        out="deepTools_qc/logs/plotCorrelation_pearson_allelic.out",
        err="deepTools_qc/logs/plotCorrelation_pearson_allelic.err"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_pearson_allelic.benchmark"
    params: 
        plotcmd = "" if plot_format == 'None' else
            "--plotFile " + "deepTools_qc/plotCorrelation/correlation.pearson.bed_coverage_allelic.heatmap." + plot_format,
        title='genes'
    shell: plotCorr_cmd


# Spearman: heatmap, scatterplot and correlation matrix
rule plotCorr_bed_spearman_allelic:
    input:
        "deepTools_qc/multiBigwigSummary/coverage_allelic.bed.npz"
    output:
        "deepTools_qc/plotCorrelation/correlation.spearman.bed_coverage_allelic.tsv"
    conda:
        CONDA_SHARED_ENV
    log:
        out="deepTools_qc/logs/plotCorrelation_spearman_allelic.out",
        err="deepTools_qc/logs/plotCorrelation_spearman_allelic.err"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_spearman_allelic.benchmark"
    params: 
        plotcmd = "" if plot_format == 'None' else
            "--plotFile " + "deepTools_qc/plotCorrelation/correlation.spearman.bed_coverage_allelic.heatmap." + plot_format,
        title='genes'
    shell: plotCorrSP_cmd


### deepTools plotPCA ##########################################################
rule plotPCA_allelic:
    input:
        "deepTools_qc/multiBigwigSummary/coverage_allelic.bed.npz"
    output:
        "deepTools_qc/plotPCA/PCA.bed_coverage_allelic.tsv"
    conda:
        CONDA_SHARED_ENV
    log:
        out="deepTools_qc/logs/plotPCA_allelic.out",
        err="deepTools_qc/logs/plotPCA_allelic.err"
    benchmark:
        "deepTools_qc/.benchmark/plotPCA_allelic.benchmark"
    params: 
        plotcmd = "" if plot_format == 'None' else
                "--plotFile " + "deepTools_qc/plotPCA/PCA.bed_coverage_allelic." + plot_format,
        title='genes'
    shell: plotPCA_cmd
