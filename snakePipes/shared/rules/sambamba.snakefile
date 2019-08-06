##sambamba is used for marking up duplications


## see Bowtie2.snakefile or RNA_mapping.snakefile for input
## takes the input from RNA mapping or DNA mapping snakefile
# mark dups
rule sambamba_markdup:
       input:
           aligner+"/{sample}.sorted.bam"
       output:
           aligner+"/{sample}.bam"# duplicate marked
       threads: 10
       log:
           out=aligner + "/logs/{sample}.sambamba_markdup.out",
           err=aligner + "/logs/{sample}.sambamba_markdup.err"
       benchmark: aligner + "/.benchmark/sambamba_markdup.{sample}.benchmark"
       conda: CONDA_SAMBAMBA_ENV
       shell: """
           sambamba markdup -t {threads} --sort-buffer-size=6000 --overflow-list-size 600000 {input} {output} 2> {log.err} > {log.out}
           """
## get statistics
rule sambamba_flagstat_sorted:
       input:
           aligner+"/{sample}.sorted.bam"
       output:
           "Sambamba/{sample}.sorted.markdup.txt"
       conda: CONDA_SAMBAMBA_ENV
       shell: """
           sambamba flagstat -p {input} > {output}
           """

rule sambamba_flagstat:
       input:
           aligner+"/{sample}.bam"
       output:
           "Sambamba/{sample}.markdup.txt"
       conda: CONDA_SAMBAMBA_ENV
       shell: """
           sambamba flagstat -p {input} > {output}
           """

## index the duplicate marked folder
rule samtools_index:
    input:
        aligner+"/{sample}.bam"
    output:
        aligner+"/{sample}.bam.bai"
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"
