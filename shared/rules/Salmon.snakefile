## Salmon Index
rule SalmonIndex:
    input:
        "Annotation/genes.filtered.fa"  ## transcripts
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
            quant = "Salmon/{sample}/quant.sf",
            quant_genes = "Salmon/{sample}/quant.genes.sf"
        benchmark:
            "Salmon/.benchmark/SalmonQuant.{sample}.benchmark"
        params:
            outdir = "Salmon/{sample}",
            libtype = salmon_libtype,
            gtf = genes_gtf,
        threads: 8
        shell:
            salmon_path+"salmon quant "
            "-p {threads} "
            "--numBootstraps 50 "
            "-g {params.gtf} "
            "-i Salmon/SalmonIndex "
            "-l {params.libtype} "
            "-1 {input.r1} -2 {input.r2} "
            "-o {params.outdir} "
else:
    rule SalmonQuant:
        input:
            fastq = fastq_dir+"/{sample}.fastq.gz",
            bin = "Salmon/SalmonIndex/sa.bin"
        output:
            quant = "Salmon/{sample}/quant.sf",
            quant_genes = "Salmon/{sample}/quant.genes.sf"
        benchmark:
            "Salmon/.benchmark/SalmonQuant.{sample}.benchmark"
        params:
            outdir = "Salmon/{sample}",
            libtype = salmon_libtype,
            gtf = genes_gtf,
        threads: 8
        shell:
            salmon_path+"salmon quant "
            "-p {threads} "
            "--numBootstraps 50 "
            "-g {params.gtf} "
            "-i Salmon/SalmonIndex "
            "-l {params.libtype} "
            "-r {input.fastq} "
            "-o {params.outdir} "


rule Salmon_symlinks:
    input:
        quant = "Salmon/{sample}/quant.sf",
        quant_genes = "Salmon/{sample}/quant.genes.sf"
    output:
        quant = "Salmon/{sample}.quant.sf",
        quant_genes = "Salmon/{sample}.quant.genes.sf"
    params:
        quant = "{sample}/quant.sf",
        quant_genes = "{sample}/quant.genes.sf"
    shell:
        "ln -s -f {params.quant} {output.quant} && touch -h {output.quant}; "
        "ln -s -f {params.quant_genes} {output.quant_genes} && touch -h {output.quant_genes} "


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


rule Salmon_genes_TPM:
    input:
        expand("Salmon/{sample}.quant.genes.sf", sample=samples)
    output:
        "Salmon/TPM.genes.tsv"
    benchmark:
        "Salmon/.benchmark/Salmon_genes_TPM.benchmark"
    log:
        "Salmon/Salmon_genes_TPM.log"
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


rule Salmon_genes_counts:
    input:
        expand("Salmon/{sample}.quant.genes.sf", sample=samples)
    output:
        "Salmon/counts.genes.tsv"
    benchmark:
        "Salmon/.benchmark/Salmon_genes_counts.benchmark"
    log:
        "Salmon/Salmon_genes_counts.log"
    shell:
        R_path+"Rscript "+os.path.join(maindir, "shared", "tools", "merge_count_tables.R")+" Name TPM {output} {input} "


## Prepare Salmon output for Sleuth
rule Salmon_wasabi:
    input:
        "Salmon/{sample}/quant.sf"
    output:
        "Salmon/{sample}/abundance.h5"
    params:
        "Salmon/{sample}"
    benchmark:
        "Salmon/.benchmark/Salmon_wasabi.benchmark"
    shell:
        "export R_LIBS_USER="+R_libs_path+" && "
        "cat "+os.path.join(workflow_tools,"wasabi.R")+" | "
        ""+os.path.join(R_path,"R")+" --vanilla --args "
        "{params}; "
