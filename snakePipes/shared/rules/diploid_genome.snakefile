from os.path import join, dirname
import glob

rule diploid_genome:
    input:
        genome = genome_fasta,
        vcf = config["VCFfile"]
    output:
        genome = "02_bwa_index/diploid_genome.fa.gz"
    threads: 4
    params:
        hap1 = strains[0],
        hap2 = strains[1]
    conda: CONDA_MAKEPAIRS_ENV
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
    input:
        genome = "02_bwa_index/diploid_genome.fa.gz"
    output:
        index = "02_bwa_index/diploid_genome.fa.gz.bwt"
    threads: 2
    resources:
        mem_mb = 50000
    conda: CONDA_MAKEPAIRS_ENV
    shell:
        """
        bwa index {input.genome}
        """

