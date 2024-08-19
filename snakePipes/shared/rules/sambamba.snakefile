##sambamba is used for marking up duplications

## see Bowtie2.snakefile or RNA_mapping.snakefile for input
## takes the input from RNA mapping or DNA mapping snakefile
# mark dups
rule sambamba_markdup:
       input:
           aligner+"/{sample}.sorted.bam"
       output:
           aligner+"/{sample}.bam"# duplicate marked
       threads: lambda wildcards: 10 if 10<max_thread else max_thread
       benchmark: aligner + "/.benchmark/sambamba_markdup.{sample}.benchmark"
       params:
           tempDir = tempDir
       conda: CONDA_SAMBAMBA_ENV
       shell: """
           TMPDIR={params.tempDir}
           MYTEMP=$(mktemp -d "${{TMPDIR:-/tmp}}"/snakepipes.XXXXXXXXXX)
           sambamba markdup -t {threads} --sort-buffer-size=6000 --overflow-list-size 600000 --tmpdir $MYTEMP {input} {output}
           rm -rf "$MYTEMP"
           """

## get statistics
rule sambamba_flagstat_sorted:
       input:
           aligner+"/{sample}.sorted.bam"
       output:
           "Sambamba/{sample}.sorted.markdup.txt"
       conda: CONDA_SAMBAMBA_ENV
       threads: lambda wildcards: 10 if 10<max_thread else max_thread
       shell: """
           sambamba flagstat -p {input} -t {threads} > {output}
           """

rule sambamba_flagstat:
       input:
           aligner+"/{sample}.bam"
       output:
           "Sambamba/{sample}.markdup.txt"
       conda: CONDA_SAMBAMBA_ENV
       threads: lambda wildcards: 10 if 10<max_thread else max_thread
       shell: """
           sambamba flagstat -p {input} -t {threads} > {output}
           """

## index the duplicate marked folder
rule samtools_index:
    input:
        aligner+"/{sample}.bam"
    output:
        aligner+"/{sample}.bam.bai"
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"

