
## featurecounts (paired options are inserted conditionally)

rule featureCounts_allele:
    input:
        gtf = "Annotation/genes.filtered.gtf",
        bam = "allelic_bams/{sample}.allele_flagged.sorted.bam",
        allele1 = "allelic_bams/{sample}.genome1.sorted.bam",
        allele2 = "allelic_bams/{sample}.genome2.sorted.bam"
    output:
        'featureCounts/{sample}.allelic_counts.txt'
    params:
        libtype = config['libraryType'],
        paired_opt = lambda wildcards: "-p -B " if config['pairedEnd'] else "",
        opts = config["featureCountsOptions"],
        tempDir = config['tempDir']
    log:
        out = "featureCounts/logs/{sample}.out",
        err = "featureCounts/logs/{sample}.err"
    threads: 8
    conda: config['CONDA_RNASEQ_ENV']
    shell: """
        TMPDIR={params.tempDir}
        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
        featureCounts \
        {params.paired_opt}{params.opts} \
        -T {threads} \
        -s {params.libtype} \
        -a {input.gtf} \
        -o {output} \
        --tmpDir $MYTEMP \
        {input.bam} {input.allele1} {input.allele2} > {log.out} 2> {log.err};
        rm -rf $MYTEMP
        """

rule merge_featureCounts:
    input:
        expand("featureCounts/{sample}.allelic_counts.txt", sample=config['samples'])
    output:
        "featureCounts/counts_allelic.tsv"
    params:
	scriptdir = os.path.join(config['maindir'], "shared", "rscripts", "merge_featureCounts.R")
    log: "featureCounts/logs/merge_featureCounts.log"
    conda: config['CONDA_RNASEQ_ENV']
    shell:
	"Rscript "+" {params.scriptdir}" +" {output} {input} 2> {log}"        
