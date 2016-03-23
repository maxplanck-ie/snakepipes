rule FASTQ_symlink:
    input:  indir+"/{sample}{read}.fastq.gz"
    output: "FASTQ/{sample}{read}.fastq.gz"
    benchmark:  "FASTQ/.benchmark/FASTQ_symlink.{sample}{read}.benchmark"
    shell:  "( [ -f {output} ] || ln -s -r {input} {output} ) && touch -h {output}"
