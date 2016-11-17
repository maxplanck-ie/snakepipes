rule bamCoverage_RPKM:
    input:
        bam = "HISAT2/{sample}.bam",
        bai = "HISAT2/{sample}.bam.bai"
    output:
        "bamCoverage/{sample}.RPKM.bw"
    params:
        bw_binsize = config["bw_binsize"]
    log:
        "bamCoverage/logs/bamCoverage_RPKM.{sample}.log"
    benchmark:
        "bamCoverage/.benchmark/bamCoverage_RPKM.{sample}.benchmark"
    threads: 8
    shell:
        deepTools_path+"bamCoverage "
        "-b {input.bam} "
        "-o {output} "
        "--binSize {params.bw_binsize} "
        "-p {threads} "
        " --normalizeUsingRPKM "
        "&> {log}"
