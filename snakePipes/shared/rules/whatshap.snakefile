rule whatshap_haplotag:
        input:
            ref = genome_fasta,
            pvcf = pvcf,
            bam = "filtered_bam/{sample}.filtered.bam",
            bai = "filtered_bam/{sample}.filtered.bam.bai"
        output:
            hbam = "allelic_bams/{sample}.allele_flagged.sorted.bam",
            hlist = "allelic_bams/{sample}_haplotype_list.tsv"
        benchmark:
            "allelic_bams/.benchmark/whatshap_haplotag.{sample}.benchmark"
        threads: 4
        conda: CONDA_WHATSHAP_ENV
        shell: """
            whatshap haplotag --ignore-read-groups -o {output.hbam} --reference {input.ref} --output-threads={threads} --output-haplotag-list={output.hlist} {input.pvcf} {input.bam}
            """

rule whatshap_split:
        input:
            hbam = "allelic_bams/{sample}.allele_flagged.sorted.bam",
            hlist = "allelic_bams/{sample}_haplotype_list.tsv"
        output:
            h1bam = "allelic_bams/{sample}.genome1.sorted.bam",
            h2bam = "allelic_bams/{sample}.genome2.sorted.bam",
            unbam = "allelic_bams/{sample}.unassigned.sorted.bam"
        benchmark:
            "allelic_bams/.benchmark/whatshap_split.{sample}.benchmark"
        conda: CONDA_WHATSHAP_ENV
        shell: """
            whatshap split  --output-h1 {output.h1bam} --output-h2 {output.h2bam} --output-untagged {output.unbam} {input.hbam} {input.hlist}
            """


rule BAMindex_allelic:
    input:
        expand("allelic_bams/{sample}.{suffix}.sorted.bam",sample=samples,suffix=['allele_flagged', 'genome1', 'genome2', 'unassigned'])
    output:
        expand("allelic_bams/{sample}.{suffix}.sorted.bam.bai",sample=samples,suffix=['allele_flagged', 'genome1', 'genome2', 'unassigned'])
    conda: CONDA_SHARED_ENV
    shell: "samtools index -M {input}"
