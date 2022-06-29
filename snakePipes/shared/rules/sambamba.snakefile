##sambamba is used for marking up duplications

## see Bowtie2.snakefile or RNA_mapping.snakefile for input
## takes the input from RNA mapping or DNA mapping snakefile
# mark dups
rule sambamba_markdup:
       input:
           config['aligner']+"/{sample}.sorted.bam"
       output:
           config['aligner']+"/{sample}.bam"# duplicate marked
       threads: lambda wildcards: 10 if 10<config['max_thread'] else config['max_thread']
       log:
           out=config['aligner'] + "/logs/{sample}.sambamba_markdup.out",
           err=config['aligner'] + "/logs/{sample}.sambamba_markdup.err"
       benchmark: config['aligner'] + "/.benchmark/sambamba_markdup.{sample}.benchmark"
       params:
           tempDir = config['tempDir']
       conda: config['CONDA_SAMBAMBA_ENV']
       shell: """
           TMPDIR={params.tempDir}
           MYTEMP=$(mktemp -d "${{TMPDIR:-/tmp}}"/snakepipes.XXXXXXXXXX)
           sambamba markdup -t {threads} --sort-buffer-size=6000 --overflow-list-size 600000 --tmpdir $MYTEMP {input} {output} 2> {log.err} > {log.out}
           rm -rf "$MYTEMP"
           """

## get statistics
rule sambamba_flagstat_sorted:
       input:
           config['aligner']+"/{sample}.sorted.bam"
       output:
           "Sambamba/{sample}.sorted.markdup.txt"
       log: "Sambamba/logs/{sample}.flagstat_sorted.log"
       conda: config['CONDA_SAMBAMBA_ENV']
       shell: """
           sambamba flagstat -p {input} > {output} 2> {log}
           """

rule sambamba_flagstat:
       input:
           config['aligner']+"/{sample}.bam"
       output:
           "Sambamba/{sample}.markdup.txt"
       log: "Sambamba/logs/{sample}.flagstat.log"
       conda: config['CONDA_SAMBAMBA_ENV']
       shell: """
           sambamba flagstat -p {input} > {output} 2> {log}
           """

## index the duplicate marked folder
rule samtools_index:
    input:
        config['aligner']+"/{sample}.bam"
    output:
        config['aligner']+"/{sample}.bam.bai"
    log: config['aligner'] + "/logs/{sample}.index.log"
    conda: config['CONDA_SHARED_ENV']
    shell: "samtools index {input} 2> {log}"
