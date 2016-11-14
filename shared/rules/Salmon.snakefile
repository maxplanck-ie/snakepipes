rule transcripts_bed2fasta:
    input:
        bed = "Annotation/genes.filtered.bed",
        genome_fasta = genome_fasta
    output:
        "Salmon/transcripts.fa"
    benchmark:
        "Salmon/.benchmark/Salmon.transcripts_bed2fasta.benchmark"
    threads: 1
    shell:
        bedtools_path+"bedtools getfasta -fi {input.genome_fasta} -bed {input.bed} -fo {output} "


##./bin/salmon index -t transcripts.fa -i transcripts_index --type quasi -k 31
