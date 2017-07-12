### deepTools bamCoverage ######################################################

rule bamCoverage:
    input:
        bam = mapping_prg+"/{sample}.bam",
        bai = mapping_prg+"/{sample}.bam.bai"
    output:
        "bamCoverage/{sample}.seq_depth_norm.bw"
    params:
        bw_binsize = bw_binsize,
        genome_size = int(genome_size),
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length)
    log:
        "bamCoverage/logs/bamCoverage.{sample}.log"
    benchmark:
        "bamCoverage/.benchmark/bamCoverage.{sample}.benchmark"
    threads: 16
    shell:
        deepTools_path+"bamCoverage "
        "-b {input.bam} "
        "-o {output} "
        "--binSize {params.bw_binsize} "
        "-p {threads} "
        "--normalizeTo1x {params.genome_size} "
        "{params.read_extension} "
        "&> {log}"


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
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
    log:
        "bamCoverage/logs/bamCoverage.{sample}.filtered.log"
    benchmark:
        "bamCoverage/.benchmark/bamCoverage.{sample}.filtered.benchmark"
    threads: 16
    shell:
        deepTools_path+"bamCoverage "
        "-b {input.bam} "
        "-o {output} "
        "--binSize {params.bw_binsize} "
        "-p {threads} "
        "--normalizeTo1x {params.genome_size} "
        "{params.read_extension} "
        "{params.blacklist} "
        "&> {log}"

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
    params:
        paired = paired,
        fragment_length = fragment_length,
        genome_size = int(genome_size),
        genome_2bit = genome_2bit,
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
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
        shell(
            deepTools_path+"computeGCBias "
            "-b {input.bam} "
            "--biasPlot {output.png} "
            "--GCbiasFrequenciesFile {output.tsv} "
            "--effectiveGenomeSize {params.genome_size} "
            "--genome {params.genome_2bit} "
            "--fragmentLength "+str(median_fragment_length)+" "
            "--sampleSize 10000000 " # very long runtime with default sample size
            "{params.blacklist} "
            "-p {threads} "
            "&> {log}"
        )


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
    shell:
        deepTools_path+"multiBamSummary bins "
        "-b {input.bams} "
        "-o {output} "
        "--labels {params.labels} "
        "--binSize 1000 "
        "{params.blacklist} "
        "-p {threads} "
        "{params.read_extension} "
        "&> {log}"


### deepTools plotCorrelation ##################################################

# Pearson: heatmap, scatterplot and correlation matrix
rule plotCorrelation_pearson:
    input:
        "deepTools_qc/multiBamSummary/read_coverage.bins.npz"
    output:
        heatpng = "deepTools_qc/plotCorrelation/correlation.pearson.read_coverage.heatmap.png",
        scatterpng = "deepTools_qc/plotCorrelation/correlation.pearson.read_coverage.scatterplot.png",
        tsv = "deepTools_qc/plotCorrelation/correlation.pearson.read_coverage.tsv"
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
        "--plotTitle 'Pearson correlation of fragment coverage' "
        "--outFileCorMatrix {output.tsv} "
        "--colorMap coolwarm "
        "--plotNumbers "
        "&> {log} && "
        +deepTools_path+"plotCorrelation "
        "-in {input} "
        "-o {output.scatterpng} "
        "--corMethod pearson "
        "--whatToPlot scatterplot "
        "--plotTitle 'Pearson correlation of fragment coverage' "
        "&>> {log}"

# Spearman: heatmap, scatterplot and correlation matrix
rule plotCorrelation_spearman:
    input:
        "deepTools_qc/multiBamSummary/read_coverage.bins.npz"
    output:
        heatpng = "deepTools_qc/plotCorrelation/correlation.spearman.read_coverage.heatmap.png",
        scatterpng = "deepTools_qc/plotCorrelation/correlation.spearman.read_coverage.scatterplot.png",
        tsv = "deepTools_qc/plotCorrelation/correlation.spearman.read_coverage.tsv"
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
        "--plotTitle 'Spearman correlation of fragment coverage' "
        "--outFileCorMatrix {output.tsv} "
        "--colorMap coolwarm "
        "--plotNumbers "
        "&> {log} && "
        +deepTools_path+"plotCorrelation "
        "-in {input} "
        "-o {output.scatterpng} "
        "--corMethod spearman "
        "--whatToPlot scatterplot "
        "--plotTitle 'Spearman correlation of fragment coverage' "
        "&>> {log}"


### deepTools plotCoverage #####################################################

rule plotCoverage:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample=samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample=samples)
    output:
        "deepTools_qc/plotCoverage/read_coverage.png"
    params:
        labels = " ".join(samples),
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length)
    log:
        "deepTools_qc/logs/plotCoverage.log"
    benchmark:
        "deepTools_qc/.benchmark/plotCoverage.benchmark"
    threads: 24
    shell:
        deepTools_path+"plotCoverage "
        "-b {input.bams} "
        "-o {output} "
        "--labels {params.labels} "
        "--plotTitle 'Genome fragment coverage without duplicates' "
        "-p {threads} "
        "{params.read_extension} "
        "--ignoreDuplicates "
        "&> {log}"


### deepTools plotPCA ##########################################################
rule plotPCA:
    input:
        "deepTools_qc/multiBamSummary/read_coverage.bins.npz"
    output:
        "deepTools_qc/plotPCA/PCA.read_coverage.png"
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
