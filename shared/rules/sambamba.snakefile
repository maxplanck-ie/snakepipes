##sambamba is used for marking up duplications

rule MarkDup:
       input:
           mapping_prg+"/{sample}.sorted.bam"
       output:
           mapping_prg+"/{sample}.bam"
       shell:
           sambamba_path+"sambamba_v0.6.6 "
           "markdup "
           "{input} "
           "{output}"

rule Flagstat:
       input:
           mapping_prg+"/{sample}.bam"
       output:
           "Sambamba/{sample}.dup.txt"
       shell:
           sambamba_path+"sambamba_v0.6.6 flagstat -p"
           " {input}"
           " | tee {output}"

