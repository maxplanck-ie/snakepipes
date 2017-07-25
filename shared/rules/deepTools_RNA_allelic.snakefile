
rule bamCoverage_RPKM_allelic:
    input:
        bam = "allelic_bams/{sample}.{suffix}.sorted.bam",
        bai = "allelic_bams/{sample}.{suffix}.sorted.bam.bai"
    output:
        "bamCoverage/allele_specific/{sample}.{suffix}.RPKM.bw"
    params:
        bw_binsize = config["bw_binsize"]
    log:
        "bamCoverage/allele_specific/logs/bamCoverage_RPKM.{sample}.{suffix}.log"
    benchmark:
        "bamCoverage/allele_specific/.benchmark/bamCoverage_RPKM.{sample}.{suffix}.benchmark"
    threads: 8
    run:
        shell(bamcov_rpkm_cmd())

rule bamCoverage_raw_allelic:
    input:
        bam = "allelic_bams/{sample}.{suffix}.sorted.bam",
        bai = "allelic_bams/{sample}.{suffix}.sorted.bam.bai"
    output:
        "bamCoverage/allele_specific/{sample}.{suffix}.coverage.bw"
    params:
        bw_binsize = bw_binsize
    log:
        "bamCoverage/allele_specificlogs/bamCoverage_coverage.{sample}.log"
    benchmark:
        "bamCoverage/allele_specific.benchmark/bamCoverage_coverage.{sample}.benchmark"
    threads: 8
    run:
        shell(bamcov_raw_cmd())

rule plotEnrichment_allelic:
    input:
        bam = expand("allelic_bams/{sample}.{suffix}.sorted.bam", sample=samples, suffix = ['genome1', 'genome2']),
        bai = expand("allelic_bams/{sample}.{suffix}.sorted.bam.bai", sample=samples, suffix = ['genome1', 'genome2']),
        bed = "Annotation/genes.filtered.bed",
    output:
        png = "deepTools_qc/plotEnrichment/plotEnrichment_allelic.png",
        tsv = "deepTools_qc/plotEnrichment/plotEnrichment_allelic.tsv",
    params:
        labels = " ".join(expand('{sample}.{suffix}', sample=samples, suffix = ['genome1', 'genome2'])),
    log:
        "deepTools_qc/logs/plotEnrichment_allelic.log"
    benchmark:
        "deepTools_qc/.benchmark/plotEnrichment_allelic.benchmark"
    threads: 8
    run:
        shell(plotEnrich_cmd())


rule multiBigwigSummary_bed_allelic:
    input:
        bw = expand("bamCoverage/allele_specific/{sample}.{suffix}.RPKM.bw", sample=samples, suffix = ['genome1', 'genome2']),
        bed = "Annotation/genes.filtered.bed",
    output:
        "deepTools_qc/multiBigwigSummary/coverage_allelic.bed.npz"
    params:
        labels = " ".join(expand('{sample}.{suffix}', sample=samples, suffix = ['genome1', 'genome2']))
    log:
        "deepTools_qc/logs/multiBigwigSummary_allelic.log"
    benchmark:
        "deepTools_qc/.benchmark/multiBigwigSummary.bed_allelic.benchmark"
    threads: 8
    run:
        shell(multiBWsum_bed_cmd())


# Pearson: heatmap, scatterplot and correlation matrix
rule plotCorr_bed_pearson_allelic:
    input:
        "deepTools_qc/multiBigwigSummary/coverage_allelic.bed.npz"
    output:
        heatpng = "deepTools_qc/plotCorrelation/correlation.pearson.bed_coverage_allelic.heatmap.png",
        scatterpng = "deepTools_qc/plotCorrelation/correlation.pearson.bed_coverage_allelic.scatterplot.png",
        tsv = "deepTools_qc/plotCorrelation/correlation.pearson.bed_coverage_allelic.tsv"
    log:
        "deepTools_qc/logs/plotCorrelation_pearson_allelic.log"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_pearson_allelic.benchmark"
    run:
        shell(plotCorr_cmd('gene'))


# Spearman: heatmap, scatterplot and correlation matrix
rule plotCorr_bed_spearman_allelic:
    input:
        "deepTools_qc/multiBigwigSummary/coverage_allelic.bed.npz"
    output:
        heatpng = "deepTools_qc/plotCorrelation/correlation.spearman.bed_coverage_allelic.heatmap.png",
        scatterpng = "deepTools_qc/plotCorrelation/correlation.spearman.bed_coverage_allelic.scatterplot.png",
        tsv = "deepTools_qc/plotCorrelation/correlation.spearman.bed_coverage_allelic.tsv"
    log:
        "deepTools_qc/logs/plotCorrelation_spearman_allelic.log"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_spearman_allelic.benchmark"
    run:
        shell(plotCorrSP_cmd('gene'))



### deepTools plotPCA ##########################################################
rule plotPCA_allelic:
    input:
        "deepTools_qc/multiBigwigSummary/coverage_allelic.bed.npz"
    output:
        "deepTools_qc/plotPCA/PCA.bed_coverage_allelic.png"
    log:
        "deepTools_qc/logs/plotPCA_allelic.log"
    benchmark:
        "deepTools_qc/.benchmark/plotPCA_allelic.benchmark"
    run:
        shell(plotPCA_cmd('gene'))
