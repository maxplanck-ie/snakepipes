rule bamCoverage_RPKM:
    input:
        bam = "HISAT2/{sample}.bam",
        bai = "HISAT2/{sample}.bam.bai"
    output:
        "BW/{sample}.RPKM.bw"
    params:
        bw_binsize = config["bw_binsize"]
    log:
        "BW/logs/bamCoverage_RPKM.{sample}.log"
    benchmark:
        "BW/.benchmark/bamCoverage_RPKM.{sample}.benchmark"
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
        bam = expand("HISAT2/{sample}.bam", sample=samples),
        bai = expand("HISAT2/{sample}.bam.bai", sample=samples),
        bed = "Annotation/genes.filtered.bed",
    output:
        png = "deepTools_qc/plotEnrichment/plotEnrichment.png",
        tsv = "deepTools_qc/plotEnrichment/plotEnrichment.tsv",
    params:
        labels = " ".join(samples),
    log:
        "deepTools_qc/plotEnrichment/plotEnrichment.log"
    benchmark:
        "deepTools_qc/.benchmark/plotEnrichment.benchmark"
    threads: 8
    shell:
        deepTools_path+"plotEnrichment "
        "-p {threads} "
        "-b {input.bam} "
        "--BED {input.bed} "
        "--plotFile {output.png} "
        "--labels {params.labels} "
        "--plotTitle 'Fraction of reads in regions' "
        "--outRawCounts {output.tsv} "
        "--variableScales "
