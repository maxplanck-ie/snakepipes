
rule bamCoverage_RPKM:
    input:
        bam = mapping_prg+"/{sample}.bam",
        bai = mapping_prg+"/{sample}.bam.bai"
    output:
        "bamCoverage/{sample}.RPKM.bw"
    params:
        bw_binsize = bw_binsize
    log:
        "bamCoverage/logs/bamCoverage_RPKM.{sample}.log"
    benchmark:
        "bamCoverage/.benchmark/bamCoverage_RPKM.{sample}.benchmark"
    threads: 8
    run:
        shell(bamcov_rpkm_cmd())

rule bamCoverage_raw:
    input:
        bam = mapping_prg+"/{sample}.bam",
        bai = mapping_prg+"/{sample}.bam.bai"
    output:
        "bamCoverage/{sample}.coverage.bw"
    params:
        bw_binsize = bw_binsize
    log:
        "bamCoverage/logs/bamCoverage_coverage.{sample}.log"
    benchmark:
        "bamCoverage/.benchmark/bamCoverage_coverage.{sample}.benchmark"
    threads: 8
    run:
        shell(bamcov_raw_cmd())

rule plotEnrichment:
    input:
        bam = expand(mapping_prg+"/{sample}.bam", sample=samples),
        bai = expand(mapping_prg+"/{sample}.bam.bai", sample=samples),
        gtf = "Annotation/genes.filtered.gtf",
        gtf2= "Annotation/genes.filtered.transcripts.gtf"
    output:
        png = "deepTools_qc/plotEnrichment/plotEnrichment.png",
        tsv = "deepTools_qc/plotEnrichment/plotEnrichment.tsv",
    params:
        labels = " ".join(samples),
    log:
        "deepTools_qc/logs/plotEnrichment.log"
    benchmark:
        "deepTools_qc/.benchmark/plotEnrichment.benchmark"
    threads: 8
    run:
        shell(plotEnrich_cmd())


rule multiBigwigSummary_bed:
    input:
        bw = expand("bamCoverage/{sample}.RPKM.bw", sample=samples),
        bed = "Annotation/genes.filtered.bed",
    output:
        "deepTools_qc/multiBigwigSummary/coverage.bed.npz"
    params:
        labels = " ".join(samples)
    log:
        "deepTools_qc/logs/multiBigwigSummary.log"
    benchmark:
        "deepTools_qc/.benchmark/multiBigwigSummary.bed.benchmark"
    threads: 8
    run:
        shell(multiBWsum_bed_cmd())


# Pearson: heatmap, scatterplot and correlation matrix
rule plotCorr_bed_pearson:
    input:
        "deepTools_qc/multiBigwigSummary/coverage.bed.npz"
    output:
        heatpng = "deepTools_qc/plotCorrelation/correlation.pearson.bed_coverage.heatmap.png",
        #scatterpng = "deepTools_qc/plotCorrelation/correlation.pearson.bed_coverage.scatterplot.png",
        tsv = "deepTools_qc/plotCorrelation/correlation.pearson.bed_coverage.tsv"
    log:
        "deepTools_qc/logs/plotCorrelation_pearson.log"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_pearson.benchmark"
    run:
        shell(plotCorr_cmd('gene'))


# Spearman: heatmap, scatterplot and correlation matrix
rule plotCorr_bed_spearman:
    input:
        "deepTools_qc/multiBigwigSummary/coverage.bed.npz"
    output:
        heatpng = "deepTools_qc/plotCorrelation/correlation.spearman.bed_coverage.heatmap.png",
        #scatterpng = "deepTools_qc/plotCorrelation/correlation.spearman.bed_coverage.scatterplot.png",
        tsv = "deepTools_qc/plotCorrelation/correlation.spearman.bed_coverage.tsv"
    log:
        "deepTools_qc/logs/plotCorrelation_spearman.log"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_spearman.benchmark"
    run:
        shell(plotCorrSP_cmd('gene'))



### deepTools plotPCA ##########################################################
rule plotPCA:
    input:
        "deepTools_qc/multiBigwigSummary/coverage.bed.npz"
    output:
        "deepTools_qc/plotPCA/PCA.bed_coverage.png"
    log:
        "deepTools_qc/logs/plotPCA.log"
    benchmark:
        "deepTools_qc/.benchmark/plotPCA.benchmark"
    run:
        shell(plotPCA_cmd('gene'))
