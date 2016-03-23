rule noDuplicates:
    input:  "Bowtie2/{sample}.bam"
    output: "noDuplicates/{sample}.noDuplicates.bam"
    benchmark:  "noDuplicates/.benchmark/noDuplicates.{sample}.benchmark"
    shell:  samtools_path+" view -hb -F1024 {input} > {output}"


rule noDuplicates_bam_index:
    input:  "noDuplicates/{sample}.noDuplicates.bam"
    output: "noDuplicates/{sample}.noDuplicates.bam.bai"
    shell:  samtools_path+" index {input}"
