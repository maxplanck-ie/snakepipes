### functions shared across workflows ##########################################
################################################################################
import os

def bamcov_cmd():
    return( (deepTools_path+"bamCoverage " +
                "-b {input.bam} " +
                "-o {output} " +
                "--binSize {params.bw_binsize} " +
                "-p {threads} " +
                "--normalizeTo1x {params.genome_size} " +
                "{params.read_extension} " +
                "&> {log}") )

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

def multiBamSum_cmd():
    return( (deepTools_path+"multiBamSummary bins " +
                    "-b {input.bams} " +
                    "-o {output} " +
                    "--labels {params.labels} " +
                    "--binSize 1000 " +
                    "{params.blacklist} " +
                    "-p {threads} " +
                    "{params.read_extension} " +
                    "&> {log}") )

def plotCorr_cmd():
    return( (deepTools_path+"plotCorrelation " +
                "-in {input} " +
                "-o {output.heatpng} " +
                "--corMethod pearson " +
                "--whatToPlot heatmap " +
                "--skipZeros " +
                "--plotTitle 'Pearson correlation of fragment coverage' " +
                "--outFileCorMatrix {output.tsv} " +
                "--colorMap coolwarm " +
                "--plotNumbers " +
                "&> {log} && " +
                deepTools_path+"plotCorrelation " +
                "-in {input} " +
                "-o {output.scatterpng} " +
                "--corMethod pearson " +
                "--whatToPlot scatterplot " +
                "--plotTitle 'Pearson correlation of fragment coverage' " +
                "&>> {log}") )


def plotCorrSP_cmd():
    return((deepTools_path+"plotCorrelation " +
        "-in {input} " +
        "-o {output.heatpng} " +
        "--corMethod spearman " +
        "--whatToPlot heatmap " +
        "--skipZeros " +
        "--plotTitle 'Spearman correlation of fragment coverage' " +
        "--outFileCorMatrix {output.tsv} " +
        "--colorMap coolwarm " +
        "--plotNumbers " +
        "&> {log} && " +
        deepTools_path+"plotCorrelation " +
        "-in {input} " +
        "-o {output.scatterpng} " +
        "--corMethod spearman " +
        "--whatToPlot scatterplot " +
        "--plotTitle 'Spearman correlation of fragment coverage' " +
        "&>> {log}"))

def plotPCA_cmd():
    return( (deepTools_path+"plotPCA " +
            "-in {input} " +
            "-o {output} " +
            "-T 'PCA of fragment coverage' " +
            "&> {log}") )

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
