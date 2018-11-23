##sambamba is used for marking up duplications


## see Bowtie2.snakefile or RNA_mapping.snakefile for input
## takes the input from RNA mapping or DNA mapping snakefile
# mark dups
rule sambamba_markdup:
       input:
           mapping_prg+"/{sample}.sorted.bam"
       output:
           mapping_prg+"/{sample}.bam"# duplicate marked
       threads: 10
       log:
           out=mapping_prg + "/logs/{sample}.sambamba_markdup.out",
           err=mapping_prg + "/logs/{sample}.sambamba_markdup.err"
       benchmark: mapping_prg + "/.benchmark/sambamba_markdup.{sample}.benchmark"
       conda: CONDA_SHARED_ENV
       shell: """
           sambamba markdup -t {threads} --sort-buffer-size=6000 --overflow-list-size 600000 {input} {output} 2> {log.err} > {log.out}
           """
## get statistics
rule sambamba_flagstat_sorted:
       input:
           mapping_prg+"/{sample}.sorted.bam"
       output:
           "Sambamba/{sample}.sorted.markdup.txt"
       conda: CONDA_SHARED_ENV
       shell: """
           sambamba flagstat -p {input} > {output}
           """

rule sambamba_flagstat:
       input:
           mapping_prg+"/{sample}.bam"
       output:
           "Sambamba/{sample}.markdup.txt"
       conda: CONDA_SHARED_ENV
       shell: """
           sambamba flagstat -p {input} > {output}
           """

## index the duplicate marked folder
rule samtools_index:
    input:
        mapping_prg+"/{sample}.bam"
    output:
        mapping_prg+"/{sample}.bam.bai"
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"
