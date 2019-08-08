
## get input bam depending on the mapping prog (use filtered bam in case of chip-seq data)
if aligner == "Bowtie2":
    rule snp_split:
        input:
            snp = SNPFile,
            bam = "filtered_bam/{sample}.filtered.bam"
        output:
            temp("filtered_bam/{sample}.filtered.sortedByName.bam"),
            expand("allelic_bams/{{sample}}.filtered.{suffix}.bam", suffix = ['allele_flagged', 'genome1', 'genome2', 'unassigned'])
        params:
            pairedEnd = '--paired' if pairedEnd else '',
            outdir = "allelic_bams"
        conda: CONDA_SHARED_ENV
        shell:
            "SNPsplit {params.pairedEnd}"
            " -o {params.outdir} --snp_file {input.snp} {input.bam}"
elif aligner == "STAR":
    rule snp_split:
        input:
            snp = SNPFile,
            bam = aligner+"/{sample}.bam"
        output:
            temp(aligner+"/{sample}.sortedByName.bam"),
            expand("allelic_bams/{{sample}}.{suffix}.bam", suffix = ['allele_flagged', 'genome1', 'genome2', 'unassigned'])
        params:
            pairedEnd = '--paired' if pairedEnd else '',
            outdir = "allelic_bams"
        conda: CONDA_SHARED_ENV
        shell:
            "SNPsplit {params.pairedEnd}"
            " -o {params.outdir} --snp_file {input.snp} {input.bam}"

# move the allele-specific bams to another folder
if aligner == "Bowtie2":
    rule movebams:
        input:
            "allelic_bams/{sample}.filtered.{suffix}.bam"
        output:
            temp("allelic_bams/{sample}.{suffix}.unsorted.bam")
        shell:
            "mv {input} {output}"
else:
    rule movebams:
        input:
            "allelic_bams/{sample}.{suffix}.bam"
        output:
            temp("allelic_bams/{sample}.{suffix}.unsorted.bam")
        shell:
            "mv {input} {output}"

# sort them
rule BAMsort_allelic:
    input: "allelic_bams/{sample}.{suffix}.unsorted.bam"
    output:
        "allelic_bams/{sample}.{suffix}.sorted.bam"
    threads:
        12
    conda: CONDA_SHARED_ENV
    shell: """
        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
        samtools sort -@ {threads} -T $MYTEMP -O bam -o {output} {input};
        rm -rf $MYTEMP
        """

# index the sorted files
rule BAMindex_allelic:
    input:
        "allelic_bams/{sample}.{suffix}.sorted.bam"
    output:
        "allelic_bams/{sample}.{suffix}.sorted.bam.bai"
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"
