rule link_bam:
    input:
        indir+"/{sample}"+bam_ext
    output:
        mapping_prg+"/{sample}"+bam_ext
    shell:
        "( [ -f {output} ] || ln -s -r {input} {output} ) && touch -h {output}"
