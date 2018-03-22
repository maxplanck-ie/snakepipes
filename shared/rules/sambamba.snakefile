##sambamba is used for marking up duplications

rule sambamba_markdup:
       input:
           mapping_prg+"/{sample}.sorted.bam"
       output:
           mapping_prg+"/{sample}.bam"
       threads: 10       
       shell:
           sambamba_path+"sambamba_v0.6.6 "
           "markdup -t {threads} --sort-buffer-size=6000 "
           "{input} "
           "{output}"

rule sambamba_flagstat:
       input:
           mapping_prg+"/{sample}.bam"
       output:
           "Sambamba/{sample}.markdup.txt"
       shell:
           sambamba_path+"sambamba_v0.6.6 flagstat -p"
           " {input}"
           " | tee {output}"


