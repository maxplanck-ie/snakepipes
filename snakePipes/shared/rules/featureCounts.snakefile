
## featurecounts (paired options are inserted conditionally)

rule featureCounts:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        gtf = "Annotation/genes.filtered.gtf",
    output:
        "featureCounts/{sample}.counts.txt"
    params:
        libtype = config['libraryType'],
        paired_opt = lambda wildcards: "-p -B " if config['pairedEnd'] else "",
        opts = config["featureCountsOptions"],
        tempDir = config['tempDir']
    log:
        out = "featureCounts/logs/{sample}.out",
        err = "featureCounts/logs/{sample}.err"
    threads: lambda wildcards: 8 if 8<config['max_thread'] else config['max_thread']
    conda: config['CONDA_RNASEQ_ENV']
    shell: """
        TMPDIR={params.tempDir}
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
        expand("featureCounts/{sample}.counts.txt", sample=config['samples'])
    output:
        "featureCounts/counts.tsv"
    params:
        scriptdir = os.path.join(config['maindir'], "shared", "rscripts", "merge_featureCounts.R")
    log: "featureCounts/logs/merge_featureCounts.log"
    conda: config['CONDA_RNASEQ_ENV']
    shell:
        "Rscript "+" {params.scriptdir}" +" {output} {input} 2> {log}"
