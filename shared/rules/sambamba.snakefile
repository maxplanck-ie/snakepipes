##sambamba is used for marking up duplications
CONDA_SHARED_ENV = "envs/shared_environment.yaml"

## see Bowtie2.snakefile or RNA_mapping.snakefile for input
## takes the input from RNA mapping or DNA mapping snakefile
# mark dups
rule sambamba_markdup:
       input:
           mapping_prg+"/{sample}.sorted.bam"
       output:
           mapping_prg+"/{sample}.bam"# duplicate marked
       threads: 10
       conda: CONDA_SHARED_ENV
       shell: """
           sambamba markdup -t {threads} --sort-buffer-size=6000 {input} {output}
           """
## get statistics
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
