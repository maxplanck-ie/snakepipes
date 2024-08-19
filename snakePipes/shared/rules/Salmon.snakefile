## Salmon Index

#rule SalmonIndex:
#    input:
#        "Annotation/genes.filtered.fa",
#        genome_fasta
#    output:
#        "Salmon/SalmonIndex/decoys.txt",
#        temp("Salmon/SalmonIndex/seq.fa"),
#        "Salmon/SalmonIndex/seq.bin"
#    benchmark:
#        "Salmon/.benchmark/Salmon.index.benchmark"
#    params:
#        salmonIndexOptions = salmonIndexOptions
#    threads: lambda wildcards: 16 if 16<max_thread else max_thread
#    conda: CONDA_RNASEQ_ENV
#    shell: """
#        grep "^>" {input[1]} | cut -d " " -f 1 | tr -d ">" > {output[0]}
#        cat {input[0]} {input[1]} > {output[1]}
#        salmon index -p {threads} -t {output[1]} -d {output[0]} -i Salmon/SalmonIndex {params.salmonIndexOptions}
#        """


def getSalmon_libtype(pairedEnd, libraryType):
    """
    Convert from a featureCounts library type to a HISAT2 option string
    """
    if pairedEnd:
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


## Salmon quant, the bootstraps are needed in Sleuth
if pairedEnd:
    rule SalmonQuant:
        input:
            r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
            r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz"
        output:
            quant = "Salmon/{sample}/quant.sf"
        benchmark:
            "Salmon/.benchmark/SalmonQuant.{sample}.benchmark"
        params:
            index = salmon_index,
            outdir = "Salmon/{sample}",
            gtf = genes_gtf,
            lib_type = getSalmon_libtype(pairedEnd, libraryType)
        threads: 8
        conda: CONDA_RNASEQ_ENV
        shell: """
            salmon quant -p {threads} --softclipOverhangs --validateMappings --numBootstraps 50 -i {params.index} -l {params.lib_type} -1 {input.r1} -2 {input.r2} -o {params.outdir}
            """
else:
    rule SalmonQuant:
        input:
            fastq = fastq_dir+"/{sample}"+reads[0]+".fastq.gz"
        output:
            quant = "Salmon/{sample}/quant.sf"
        benchmark:
            "Salmon/.benchmark/SalmonQuant.{sample}.benchmark"
        params:
            index = salmon_index,
            outdir = "Salmon/{sample}",
            gtf = genes_gtf,
            lib_type = getSalmon_libtype(pairedEnd, libraryType)
        threads: 8
        conda: CONDA_SALMON_ENV
        shell: """
            salmon quant -p {threads} --softclipOverhangs --validateMappings --numBootstraps 50 -i {params.index} -l {params.lib_type} -r {input.fastq} -o {params.outdir}
            """


rule Salmon_symlinks:
    input:
        quant = "Salmon/{sample}/quant.sf"
    output:
        quant = "Salmon/{sample}.quant.sf"
    shell: """
         ln -s ../{input.quant} {output.quant}
            """


rule Salmon_TPM:
    input:
        expand("Salmon/{sample}.quant.sf", sample=samples)
    output:
        "Salmon/TPM.transcripts.tsv"
    benchmark:
        "Salmon/.benchmark/Salmon_TPM.benchmark"
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
