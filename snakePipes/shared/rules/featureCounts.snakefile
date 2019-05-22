
## featurecounts (paired options are inserted conditionally)

rule featureCounts:
    input:
        bam = mapping_prg+"/{sample}.bam",
        gtf = "Annotation/genes.filtered.gtf",
    output:
        "featureCounts/{sample}.counts.txt"
    params:
        libtype = library_type,
        paired_opt = lambda wildcards: "-p -B " if paired else "",
        opts = config["featurecounts_options"],
    log:
        out = "featureCounts/{sample}.out",
        err = "featureCounts/{sample}.err"
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
    conda: CONDA_RNASEQ_ENV
    shell:
        "Rscript "+os.path.join(maindir, "shared", "rscripts", "merge_featureCounts.R")+" {output} {input}"
