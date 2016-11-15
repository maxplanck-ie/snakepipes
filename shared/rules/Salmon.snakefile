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


rule SalmonIndex:
    input:
        "Salmon/transcripts.fa"
    output:
        "Salmon/SalmonIndex/sa.bin"
    benchmark:
        "Salmon/.benchmark/Salmon.index.benchmark"
    params:
        salmon_index_options = salmon_index_options
    log: "Salmon/SalmonIndex/SalmonIndex.log"
    threads: 8
    shell:
        salmon_path+" index -p {threads} -t {input} -i Salmon/SalmonIndex {params.salmon_index_options} &> {log} && touch {output}"

## Salmon quant
if paired:
    rule SalmonQuant:
        input:
            r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
            r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz",
            bin = "Salmon/SalmonIndex/sa.bin"
        output:
            "Salmon/transcripts_quant/Salmon.quant.{sample}.log"
        benchmark:
            "Salmon/.benchmark/Salmon.quant.{sample}.benchmark"
        params:
            outdir = "Salmon/transcripts_quant/{sample}"
        threads: 8
        shell:
            salmon_path+" quant -p {threads} -i Salmon/SalmonIndex -l a -1 {input.r1} -2 {input.r2} -o {params.outdir} &> {output}"

else:
    rule SalmonQuant:
        input:
            fastq = fastq_dir+"/{sample}.fastq.gz",
            bin = "Salmon/SalmonIndex/sa.bin"
        output:
            "Salmon/transcripts_quant/Salmon.quant.{sample}.log"
        benchmark:
            "Salmon/.benchmark/Salmon.quant.{sample}.benchmark"
        params:
            outdir = "Salmon/transcripts_quant/{sample}"
        threads: 8
        shell:
            salmon_path+" quant -p {threads} -i Salmon/SalmonIndex -l a -r {input} -o {params.outdir} &> {output}"
