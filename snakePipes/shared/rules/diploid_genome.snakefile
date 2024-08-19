from os.path import join, dirname
import glob


rule diploid_genome:
    # given a genome_fasta and a vcf-file with two haplotypes
    # create two haplotype genomes and concatenate to
    # one fasta
    input:
        genome=genome_fasta,
        vcf=config["VCFfile"],
    output:
        genome="genome/diploid_genome.fa.gz",
    threads: 4
    benchmark:
        "genome/.benchmark/diploid_genome.benchmark"
    params:
        hap1=strains[0],
        hap2=strains[1],
    conda:
        CONDA_MAKEPAIRS_ENV
    shell:
        """
        bcftools consensus \
            --fasta-ref {input.genome} \
            --haplotype 1 {input.vcf} \
            --sample {params.hap1} \
            | sed -E 's/(>[^[:space:]]+).*/\\1_{params.hap1}/g' \
            | bgzip -c > genome_{params.hap1}.fa.gz

        bcftools consensus \
            --fasta-ref {input.genome} \
            --haplotype 1 {input.vcf} \
            --sample {params.hap2} \
            | sed -E 's/(>[^[:space:]]+).*/\\1_{params.hap2}/g' \
            | bgzip -c > genome_{params.hap2}.fa.gz

        cat genome_{params.hap1}.fa.gz genome_{params.hap2}.fa.gz \
            > {output.genome}

        rm genome_{params.hap1}.fa.gz genome_{params.hap2}.fa.gz
        """


rule bwa_index_diploid_genome:
    # index concatenated genome generated by rule diploid_genome
    input:
        genome="genome/diploid_genome.fa.gz",
    output:
        index="genome/diploid_genome.fa.gz.bwt",
    threads: 2
    resources:
        mem_mb=50000,
    conda:
        CONDA_MAKEPAIRS_ENV
    shell:
        """
        bwa index {input.genome}
        """


rule chr_sizes:
    # obtain chromosome sizes from bwa index
    input:
        bwaix="genome/{ref}.fa.gz.bwt",
    output:
        chromsize="genome/{ref}.chromsizes",
    params:
        fnagz=lambda wildcards, input: Path(input.bwaix).with_suffix(""),
    threads: 1
    conda:
        CONDA_MAKEPAIRS_ENV
    shell:
        """
        samtools faidx {params.fnagz}
        cut -f1,2 {params.fnagz}.fai > {output.chromsize}
        """
