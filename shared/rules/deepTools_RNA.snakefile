rule bamCoverage_RPKM:
    input:
        bam = mapping_prg+"/{sample}.bam",
        bai = mapping_prg+"/{sample}.bam.bai"
    output:
        "Tracks/{sample}.RPKM.bw"
    params:
        bw_binsize = bw_binsize
    log:
        "Tracks/logs/bamCoverage_RPKM.{sample}.log"
    benchmark:
        "Tracks/.benchmark/bamCoverage_RPKM.{sample}.benchmark"
    threads: 8
    shell:
        deepTools_path+"bamCoverage "
        "-b {input.bam} "
        "-o {output} "
        "--binSize {params.bw_binsize} "
        "-p {threads} "
        " --normalizeUsingRPKM "
        "&> {log}"


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
    shell:
        deepTools_path+"plotEnrichment "
        "-p {threads} "
        "-b {input.bam} "
        "--BED {input.gtf} {input.gtf2} "
        "--plotFile {output.png} "
        "--labels {params.labels} "
        "--plotTitle 'Fraction of reads in regions' "
        "--outRawCounts {output.tsv} "
        "&> {log} "


rule multiBigwigSummary_bed:
    input:
        bw = expand("BW/{sample}.RPKM.bw", sample=samples),
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
    shell:
        deepTools_path+"multiBigwigSummary BED-file "
        "--BED {input.bed} "
        "-b {input.bw} "
        "-o {output} "
        "--labels {params.labels} "
        "--binSize 1000 "
        "-p {threads} "
        "&> {log} "


# Pearson: heatmap, scatterplot and correlation matrix
rule plotCorrelation_bed_pearson:
    input:
        "deepTools_qc/multiBigwigSummary/coverage.bed.npz"
    output:
        heatpng = "deepTools_qc/plotCorrelation/correlation.pearson.bed_coverage.heatmap.png",
        scatterpng = "deepTools_qc/plotCorrelation/correlation.pearson.bed_coverage.scatterplot.png",
        tsv = "deepTools_qc/plotCorrelation/correlation.pearson.bed_coverage.tsv"
    log:
        "deepTools_qc/logs/plotCorrelation_pearson.log"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_pearson.benchmark"
    shell:
        deepTools_path+"plotCorrelation "
        "-in {input} "
        "-o {output.heatpng} "
        "--corMethod pearson "
        "--whatToPlot heatmap "
        "--skipZeros "
        "--plotTitle 'Pearson correlation of region coverage' "
        "--outFileCorMatrix {output.tsv} "
        "--colorMap coolwarm "
        "--plotNumbers "
        "&> {log} && "
        +deepTools_path+"plotCorrelation "
        "-in {input} "
        "-o {output.scatterpng} "
        "--corMethod pearson "
        "--whatToPlot scatterplot "
        "--plotTitle 'Pearson correlation of region coverage' "
        "&>> {log}"


# Spearman: heatmap, scatterplot and correlation matrix
rule plotCorrelation_bed_spearman:
    input:
        "deepTools_qc/multiBigwigSummary/coverage.bed.npz"
    output:
        heatpng = "deepTools_qc/plotCorrelation/correlation.spearman.bed_coverage.heatmap.png",
        scatterpng = "deepTools_qc/plotCorrelation/correlation.spearman.bed_coverage.scatterplot.png",
        tsv = "deepTools_qc/plotCorrelation/correlation.spearman.bed_coverage.tsv"
    log:
        "deepTools_qc/logs/plotCorrelation_spearman.log"
    benchmark:
        "deepTools_qc/.benchmark/plotCorrelation_spearman.benchmark"
    shell:
        deepTools_path+"plotCorrelation "
        "-in {input} "
        "-o {output.heatpng} "
        "--corMethod spearman "
        "--whatToPlot heatmap "
        "--skipZeros "
        "--plotTitle 'Spearman correlation of region coverage' "
        "--outFileCorMatrix {output.tsv} "
        "--colorMap coolwarm "
        "--plotNumbers "
        "&> {log} && "
        +deepTools_path+"plotCorrelation "
        "-in {input} "
        "-o {output.scatterpng} "
        "--corMethod spearman "
        "--whatToPlot scatterplot "
        "--plotTitle 'Spearman correlation of region coverage' "
        "&>> {log}"



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
    shell:
        deepTools_path+"plotPCA "
            "-in {input} "
            "-o {output} "
            " -T 'PCA of fragment coverage' "
            "&> {log}"
