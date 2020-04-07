
## featurecounts (paired options are inserted conditionally)

rule featureCounts:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        gtf = "Annotation/genes.filtered.gtf",
    output:
        "featureCounts/{sample}.counts.txt"
    params:
        libtype = libraryType,
        paired_opt = lambda wildcards: "-p -B " if pairedEnd else "",
        opts = config["featureCountsOptions"],
    log:
        out = "featureCounts/logs/{sample}.out",
        err = "featureCounts/logs/{sample}.err"
    threads: 8
    conda: CONDA_RNASEQ_ENV
    shell: """
        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
        featureCounts  \
        {params.paired_opt}{params.opts} \
        -T {threads} \
        -s {params.libtype} \
        -a {input.gtf} \
        -o {output} \
        --tmpDir $MYTEMP \
        {input.bam} > {log.out} 2> {log.err};
        rm -rf $MYTEMP
        """

rule merge_featureCounts:
    input:
        expand("featureCounts/{sample}.counts.txt", sample=samples)
    output:
        "featureCounts/counts.tsv"
    log: "featureCounts/logs/merge_featureCounts.log"
    conda: CONDA_RNASEQ_ENV
    shell:
        "Rscript "+os.path.join(maindir, "shared", "rscripts", "merge_featureCounts.R")+" {output} {input} 2> {log}"
