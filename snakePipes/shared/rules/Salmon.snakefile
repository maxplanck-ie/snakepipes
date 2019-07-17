## Salmon Index
rule SalmonIndex:
    input:
        "Annotation/genes.filtered.fa"  ## transcripts
    output:
        "Salmon/SalmonIndex/sa.bin"
    benchmark:
        "Salmon/.benchmark/Salmon.index.benchmark"
    params:
        salmonIndexOptions = salmonIndexOptions
    log:
        out = "Salmon/SalmonIndex/SalmonIndex.out",
        err = "Salmon/SalmonIndex/SalmonIndex.err",
    threads: 16
    conda: CONDA_RNASEQ_ENV
    shell: """
        salmon index -p {threads} -t {input} -i Salmon/SalmonIndex {params.salmonIndexOptions} > {log.out} 2> {log.err}
        """


def getSalmon_libtype(paired, libraryType):
    """
    Convert from a featureCounts library type to a HISAT2 option string
    """
    if paired:
        if libraryType == 1:
            return "ISF"
        elif libraryType == 2:
            return "ISR"
        else:
            return "IU"
    else:
        if libraryType == 1:
            return "SF"
        elif libraryType == 2:
            return "SR"
        else:
            return "U"


## Salmon quant
if paired:
    rule SalmonQuant:
        input:
            r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
            r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz",
            bin = "Salmon/SalmonIndex/sa.bin",
        output:
            quant = "Salmon/{sample}/quant.sf",
            quant_genes = "Salmon/{sample}/quant.genes.sf"
        benchmark:
            "Salmon/.benchmark/SalmonQuant.{sample}.benchmark"
        params:
            outdir = "Salmon/{sample}",
            gtf = genes_gtf,
            lib_type = getSalmon_libtype(paired, libraryType)
        threads: 8
        conda: CONDA_RNASEQ_ENV
        shell: """
            salmon quant -p {threads} --numBootstraps 50 -g {params.gtf} -i Salmon/SalmonIndex -l {params.lib_type} -1 {input.r1} -2 {input.r2} -o {params.outdir}
            """
else:
    rule SalmonQuant:
        input:
            fastq = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
            bin = "Salmon/SalmonIndex/sa.bin",
        output:
            quant = "Salmon/{sample}/quant.sf",
            quant_genes = "Salmon/{sample}/quant.genes.sf"
        benchmark:
            "Salmon/.benchmark/SalmonQuant.{sample}.benchmark"
        params:
            outdir = "Salmon/{sample}",
            gtf = genes_gtf,
            lib_type = getSalmon_libtype(paired, libraryType)
        threads: 8
        conda: CONDA_RNASEQ_ENV
        shell: """
            salmon quant -p {threads} --numBootstraps 50 -g {params.gtf} -i Salmon/SalmonIndex -l {params.lib_type} -r {input.fastq} -o {params.outdir}
            """


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
    shell: """
        ln -s -f {params.quant} {output.quant}
        ln -s -f {params.quant_genes} {output.quant_genes}
        """


rule Salmon_TPM:
    input:
        expand("Salmon/{sample}.quant.sf", sample=samples)
    output:
        "Salmon/TPM.transcripts.tsv"
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
        "Salmon/counts.transcripts.tsv"
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
        "Rscript "+os.path.join(maindir, "shared", "rscripts", "merge_count_tables.R")+" Name NumReads {output} {input} "


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
