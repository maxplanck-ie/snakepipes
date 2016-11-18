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
        salmon_path+"salmon index -p {threads} -t {input} -i Salmon/SalmonIndex {params.salmon_index_options} &> {log} && touch {output}"


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
            gtf = genes_gtf,
        threads: 8
        shell:
            salmon_path+"salmon quant -p {threads} -g {params.gtf} -i Salmon/SalmonIndex -l {params.libtype} -1 {input.r1} -2 {input.r2} -o {params.outdir}"

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
            gtf = genes_gtf,
        threads: 8
        shell:
            salmon_path+"salmon quant -p {threads} -g {params.gtf} -i Salmon/SalmonIndex -l {params.libtype} -r {input.fastq} -o {params.outdir}"


rule Salmon_symlinks:
    input:
        expand("Salmon/{sample}/quant.sf", sample=samples)
    output:
        "Salmon/{sample}.quant.sf"
    shell:
        "( [ -f {output} ] || ln -s -r {input} {output} ) && touch -h {output}"


rule Salmon_TPM:
    input:
        expand("Salmon/{sample}.quant.sf", sample=samples)
    output:
        "Salmon/TPM.tsv"
    benchmark:
        "Salmon/.benchmark/Salmon_TPM.benchmark"
    log:
        "Salmon/Salmon_TPM.log"
    shell:
        R_path+"Rscript "+os.path.join(maindir, "shared", "tools", "merge_count_tables.R")+" Name TPM {output} {input} "


rule Salmon_counts:
    input:
        expand("Salmon/{sample}.quant.sf", sample=samples)
    output:
        "Salmon/counts.tsv"
    benchmark:
        "Salmon/.benchmark/Salmon_counts.benchmark"
    log:
        "Salmon/Salmon_counts.log"
    shell:
        R_path+"Rscript "+os.path.join(maindir, "shared", "tools", "merge_count_tables.R")+" Name NumReads {output} {input} "
