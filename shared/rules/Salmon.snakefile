## Salmon Index
rule SalmonIndex:
    input:
        "Annotation/genes_filtered.fa"  ## transcripts
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
            "Salmon/{sample}/quant.sf"
        benchmark:
            "Salmon/.benchmark/SalmonQuant.{sample}.benchmark"
        params:
            outdir = "Salmon/{sample}",
            libtype = salmon_libtype,
        threads: 8
        shell:
            salmon_path+" quant -p {threads} -i Salmon/SalmonIndex -l {params.libtype} -1 {input.r1} -2 {input.r2} -o {params.outdir}"

else:
    rule SalmonQuant:
        input:
            fastq = fastq_dir+"/{sample}.fastq.gz",
            bin = "Salmon/SalmonIndex/sa.bin"
        output:
            "Salmon/{sample}/quant.sf"
        benchmark:
            "Salmon/.benchmark/SalmonQuant.{sample}.benchmark"
        params:
            outdir = "Salmon/{sample}",
            libtype = salmon_libtype,
        threads: 8
        shell:
            salmon_path+" quant -p {threads} -i Salmon/SalmonIndex -l {params.libtype} -r {input.fastq} -o {params.outdir}"


rule Salmon_get_TPMs:
    input:
        "Salmon/{sample}/quant.sf"
    output:
        temp("Salmon/{sample}.TPM.tsv")
    shell:
        "tail -n +2 {input} | cut -f4 > {output}"


rule Salmon_merge_TPMs:
    input:
        expand("Salmon/{sample}.TPM.tsv", sample=samples)
    output:
        "Salmon/TPM.tsv"
    params:
        "\t".join(samples)
    shell:
        "echo '{params}' > {output} && paste {input} >> {output}"


rule Salmon_get_counts:
    input:
        "Salmon/{sample}/quant.sf"
    output:
        temp("Salmon/{sample}.counts.tsv")
    shell:
        "tail -n +2 {input} | cut -f5 > {output}"


rule Salmon_merge_counts:
    input:
        expand("Salmon/{sample}.counts.tsv", sample=samples)
    output:
        "Salmon/counts.tsv"
    params:
        "\t".join(samples)
    shell:
        "echo '{params}' > {output} && paste {input} >> {output}"
