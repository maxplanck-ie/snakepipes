rule link_bam:
    input:
        indir + "/{sample}" + bamExt
    output:
        mapping_prg + "/{sample}.bam"
    shell:
        "( [ -f {output} ] || ln -s -r {input} {output} )"


rule samtools_index_external:
    input:
        mapping_prg + "/{sample}.bam"
    output:
        mapping_prg + "/{sample}.bam.bai"
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"


rule link_bam_bai_external:
    input:
        bam = mapping_prg + "/{sample}.bam",
        bai = mapping_prg + "/{sample}.bam.bai"
    output:
        bam_out = "filtered_bam/{sample}.filtered.bam",
        bai_out = "filtered_bam/{sample}.filtered.bam.bai",
    shell:
        "( [ -f {output.bam_out} ] || ( ln -s -r {input.bam} {output.bam_out} && ln -s -r {input.bai} {output.bai_out} ) )"


rule sambamba_flagstat:
       input:
           mapping_prg + "/{sample}.bam"
       output:
           "Sambamba/{sample}.markdup.txt"
       conda: CONDA_SAMBAMBA_ENV
       shell: """
           sambamba flagstat -p {input} > {output}
           """
