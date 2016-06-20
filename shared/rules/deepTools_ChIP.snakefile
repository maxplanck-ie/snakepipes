### plotFingerprint ############################################################

rule plotFingerprint:
    input:
        bams = expand("Bowtie2/{sample}.bam", sample=samples),
        bais = expand("Bowtie2/{sample}.bam.bai", sample=samples)
    output:
        "deepTools_ChIP/plotFingerprint/fingerprint.png"
    params:
        labels = " ".join(samples),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed else "",
        read_extension = "--extendReads" if paired else "--extendReads "+str(fragment_length)
    log:
        "deepTools_ChIP/logs/plotFingerprint.log"
    benchmark:
        "deepTools_ChIP/.benchmark/plotFingerprint.benchmark"
    threads: 32
    shell:
        deepTools_path+"plotFingerprint "
        "-b {input.bams} "
        "--plotFile {output} "
        "--labels {params.labels} "
        "--plotTitle 'Cumulative read counts per bin' "
        "{params.blacklist} "
        "-p {threads} "
        "{params.read_extension} "
        "&> {log}"
