
rule convertLibraryTypeSalmon:
    input: "Annotation/genes.filtered.fa"
    output: "Salmon/lib_type.txt"
    params:
        lib_str = "PE" if paired else "SE",
        from_library_type = library_type,
        from_prg = "featureCounts",
        to_prg="Salmon",
        tsv = os.path.join(maindir, "shared", "tools", "library_type.tsv"),
        rscript = os.path.join(maindir, "shared", "rscripts", "library_type.R")
    threads: 1
    conda: CONDA_RNASEQ_ENV
    shell:
        "Rscript {params.rscript} {params.tsv} {params.lib_str} {params.from_library_type} {params.from_prg} {params.to_prg} > {output}"


## Salmon Index
rule SalmonIndex:
    input:
        "Annotation/genes.filtered.fa"  ## transcripts
    output:
        "Salmon/SalmonIndex/sa.bin"
    benchmark:
        "Salmon/.benchmark/Salmon.index.benchmark"
    params:
        salmon_index_options = salmon_index_options,
        sasamp = 1

    log:
        out = "Salmon/SalmonIndex/SalmonIndex.out",
        err = "Salmon/SalmonIndex/SalmonIndex.err",
    threads: 8
    conda: CONDA_RNASEQ_ENV
    shell: """
        salmon index --sasamp {params.sasamp} -p {threads} -t {input} -i Salmon/SalmonIndex {params.salmon_index_options} > {log.out} 2> {log.err}
        """
#        touch {output}

## Salmon quant
if paired:
    rule SalmonQuant:
        input:
            r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
            r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz",
            bin = "Salmon/SalmonIndex/sa.bin",
            lib_type = "Salmon/lib_type.txt"
        output:
            quant = "Salmon/{sample}/quant.sf",
            quant_genes = "Salmon/{sample}/quant.genes.sf"
        benchmark:
            "Salmon/.benchmark/SalmonQuant.{sample}.benchmark"
        params:
            outdir = "Salmon/{sample}",
            gtf = genes_gtf,
        threads: 8
        conda: CONDA_RNASEQ_ENV
        shell:
            "lib_type=$(cat {input.lib_type} ); echo \"lib_type=\"$lib_type 1>&2; "
            "salmon quant "
            "-p {threads} "
            "--numBootstraps 50 "
            "-g {params.gtf} "
            "-i Salmon/SalmonIndex "
            "-l $lib_type "
            "-1 {input.r1} -2 {input.r2} "
            "-o {params.outdir} "
else:
    rule SalmonQuant:
        input:
            fastq = fastq_dir+"/{sample}.fastq.gz",
            bin = "Salmon/SalmonIndex/sa.bin",
            lib_type = "Salmon/lib_type.txt"
        output:
            quant = "Salmon/{sample}/quant.sf",
            quant_genes = "Salmon/{sample}/quant.genes.sf"
        benchmark:
            "Salmon/.benchmark/SalmonQuant.{sample}.benchmark"
        params:
            outdir = "Salmon/{sample}",
            gtf = genes_gtf,
        threads: 8
        conda: CONDA_RNASEQ_ENV
        shell:
            "lib_type=$(cat {input.lib_type} ); echo \"lib_type=\"$lib_type 1>&2; "
            "salmon quant "
            "-p {threads} "
            "--numBootstraps 50 "
            "-g {params.gtf} "
            "-i Salmon/SalmonIndex "
            "-l $lib_type "
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
        "ln -s -f {params.quant} {output.quant}" #&& touch -h {output.quant};
        "ln -s -f {params.quant_genes} {output.quant_genes}" # && touch -h {output.quant_genes}


rule Salmon_TPM:
    input:
        expand("Salmon/{sample}.quant.sf", sample=samples)
    output:
        "Salmon/TPM.tsv"
    benchmark:
        "Salmon/.benchmark/Salmon_TPM.benchmark"
    log:
        "Salmon/Salmon_TPM.log"
    conda: CONDA_RNASEQ_ENV
    shell:
        "Rscript "+os.path.join(maindir, "shared", "rscripts", "merge_count_tables.R")+" Name TPM {output} {input} "


rule Salmon_genes_TPM:
    input:
        expand("Salmon/{sample}.quant.genes.sf", sample=samples)
    output:
        "Salmon/TPM.genes.tsv"
    benchmark:
        "Salmon/.benchmark/Salmon_genes_TPM.benchmark"
    log:
        "Salmon/Salmon_genes_TPM.log"
    conda: CONDA_RNASEQ_ENV
    shell:
        "Rscript "+os.path.join(maindir, "shared", "rscripts", "merge_count_tables.R")+" Name TPM {output} {input} "


rule Salmon_counts:
    input:
        expand("Salmon/{sample}.quant.sf", sample=samples)
    output:
        "Salmon/counts.tsv"
    benchmark:
        "Salmon/.benchmark/Salmon_counts.benchmark"
    log:
        "Salmon/Salmon_counts.log"
    conda: CONDA_RNASEQ_ENV
    shell:
        "Rscript "+os.path.join(maindir, "shared", "rscripts", "merge_count_tables.R")+" Name NumReads {output} {input} "


rule Salmon_genes_counts:
    input:
        expand("Salmon/{sample}.quant.genes.sf", sample=samples)
    output:
        "Salmon/counts.genes.tsv"
    benchmark:
        "Salmon/.benchmark/Salmon_genes_counts.benchmark"
    log:
        "Salmon/Salmon_genes_counts.log"
    conda: CONDA_RNASEQ_ENV
    shell:
        "Rscript "+os.path.join(maindir, "shared", "rscripts", "merge_count_tables.R")+" Name TPM {output} {input} "


## Prepare Salmon output for Sleuth
rule Salmon_wasabi:
    input:
        "Salmon/{sample}.quant.sf"
    output:
        "Salmon/{sample}/abundance.h5"
    params:
        "Salmon/{sample}/"
    conda: CONDA_RNASEQ_ENV
    shell:
        "Rscript "+os.path.join(maindir, "shared", "rscripts", "wasabi.R")+" {params}"
