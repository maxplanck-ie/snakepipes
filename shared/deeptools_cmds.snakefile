### functions shared across workflows ##########################################
################################################################################
import os

# bamCoverage RAW
def bamcov_raw_cmd():
    return(deepTools_path+"bamCoverage " +
            "-b {input.bam} " +
            "-o {output} " +
            "--binSize {params.bw_binsize} " +
            "-p {threads} " +
            "&> {log}")

# bamCoverage CHIP
def bamcov_cmd():
    return(deepTools_path+"bamCoverage " +
                "-b {input.bam} " +
                "-o {output} " +
                "--binSize {params.bw_binsize} " +
                "-p {threads} " +
                "--normalizeTo1x {params.genome_size} " +
                "--ignoreForNormalization {params.ignore} " +
                "{params.read_extension} " +
                "&> {log}")

## bamCoverage RNAseq
def bamcov_rpkm_cmd():
    return( (deepTools_path+"bamCoverage "
            "-b {input.bam} "
            "-o {output} "
            "--binSize {params.bw_binsize} "
            "-p {threads} "
            " --normalizeUsingRPKM "
            "&> {log}") )

## computeGC bias (chIPseq)
def gcbias_cmd():
            if params.paired:
                median_fragment_length = cf.get_fragment_length(input.insert_size_metrics)
            else:
                median_fragment_length = params.fragment_length
            shell(
                (deepTools_path+"computeGCBias " +
                "-b {input.bam} " +
                "--biasPlot {output.png} " +
                "--GCbiasFrequenciesFile {output.tsv} " +
                "--effectiveGenomeSize {params.genome_size} " +
                "--genome {params.genome_2bit} " +
                "--fragmentLength "+str(median_fragment_length)+" " +
                "--sampleSize 10000000 " # very long runtime with default sample size
                "{params.blacklist} " +
                "-p {threads} " +
                "&> {log}")
                )

# plot Enrichment (RNAseq)
def plotEnrich_cmd():
    return( (deepTools_path+"plotEnrichment "
        "-p {threads} "
        "-b {input.bam} "
        "--BED {input.gtf} {input.gtf2} "
        "--plotFile {output.png} "
        "--labels {params.labels} "
        "--plotTitle 'Fraction of reads in regions' "
        "--outRawCounts {output.tsv} "
        "--variableScales "
        "&> {log}") )

# multiBAMsum ChIP
def multiBamSummary_cmd():
    return( (deepTools_path+"multiBamSummary bins " +
                    "-b {input.bams} " +
                    "-o {output} " +
                    "--labels {params.labels} " +
                    "--binSize 1000 " +
                    "{params.blacklist} " +
                    "-p {threads} " +
                    "{params.read_extension} " +
                    "&> {log}") )

# multiBAMsum RNA
def multiBWsum_bed_cmd():
    return( (deepTools_path+"multiBigwigSummary BED-file "
                "--BED {input.bed} "
                "-b {input.bw} "
                "-o {output} "
                "--labels {params.labels} "
                "--binSize 1000 "
                "-p {threads} "
                "&> {log} "))

## plot Corr (both)
def plotCorr_cmd(what):
    return( (deepTools_path+"plotCorrelation " +
                "-in {input} " +
                "-o {output.heatpng} " +
                "--corMethod pearson " +
                "--whatToPlot heatmap " +
                "--skipZeros " +
                "--plotTitle 'Pearson correlation of "+what+" coverage' " +
                "--outFileCorMatrix {output.tsv} " +
                "--colorMap coolwarm " +
                "--plotNumbers " +
                "&> {log} && " +
                deepTools_path+"plotCorrelation " +
                "-in {input} " +
                "-o {output.scatterpng} " +
                "--corMethod pearson " +
                "--whatToPlot scatterplot " +
                "--plotTitle 'Pearson correlation of "+what+" coverage' " +
                "&>> {log}") )

## plot Corr Spearman (both)
def plotCorrSP_cmd(what):
    return( (deepTools_path+"plotCorrelation " +
        "-in {input} " +
        "-o {output.heatpng} " +
        "--corMethod spearman " +
        "--whatToPlot heatmap " +
        "--skipZeros " +
        "--plotTitle 'Spearman correlation of "+what+" coverage' " +
        "--outFileCorMatrix {output.tsv} " +
        "--colorMap coolwarm " +
        "--plotNumbers " +
        "&> {log} && " +
        deepTools_path+"plotCorrelation " +
        "-in {input} " +
        "-o {output.scatterpng} " +
        "--corMethod spearman " +
        "--whatToPlot scatterplot " +
        "--plotTitle 'Spearman correlation of "+what+" coverage' " +
        "&>> {log}"))

# plot PCA (both)
def plotPCA_cmd(what):
    return( (deepTools_path+"plotPCA " +
            "-in {input} " +
            "-o {output} " +
            "-T 'PCA of fragment coverage' " +
            "&> {log}") )

# plot Coverage
def plotCov_cmd():
    return( (deepTools_path+"plotCoverage " +
                "-b {input.bams} " +
                "-o {output} " +
                "--labels {params.labels} " +
                "--plotTitle 'Genome fragment coverage without duplicates' " +
                "-p {threads} " +
                "{params.read_extension} " +
                "--ignoreDuplicates " +
                "&> {log}") )
