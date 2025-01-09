checkpoint whatshap_haplotag:
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

checkpoint whatshap_split:
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

#wildcard_constraints:
#    suffix = ['allele_flagged', 'genome1', 'genome2', 'unassigned']


def collect_split_bams(wildcards):
      checkpoint_output_a = checkpoints.whatshap_haplotag.get(**wildcards).output["hbam"]
      checkpoint_output_b = checkpoints.whatshap_split.get(**wildcards).output
      checkpoint_output = checkpoint_output_a + checkpoint_output_b
      return expand("allelic_bams/{{sample}}.{suffix}.sorted.bam",
                  suffix = glob_wildcards("allelic_bams/{sample}.{suffix}.sorted.bam").suffix)

rule BAMindex_allelic:
    input:
#        "allelic_bams/{sample}.{suffix}.sorted.bam"
        collect_split_bams
    output:
        "allelic_bams/{sample}.{suffix}.sorted.bam.bai"
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"
