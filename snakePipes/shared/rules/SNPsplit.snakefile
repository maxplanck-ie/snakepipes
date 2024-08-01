
## get input bam depending on the mapping prog (use filtered bam in case of chipseq data)
if aligner == "Bowtie2":
    rule snp_split:
        input:
            snp = SNPFile,
            bam = "filtered_bam/{sample}.filtered.bam"
        output:
            targetbam = expand("allelic_bams/{{sample}}.filtered.{suffix}.bam", suffix = ['allele_flagged', 'genome1', 'genome2', 'unassigned']),
            tempbam = temp("filtered_bam/{sample}.filtered.sortedByName.bam"),
            rep1 = "allelic_bams/{sample}.filtered.SNPsplit_report.yaml",
            rep2 = "allelic_bams/{sample}.filtered.SNPsplit_sort.yaml"
        log: "allelic_bams/logs/{sample}.snp_split.log"
        params:
            pairedEnd = '--paired' if pairedEnd else '',
            outdir = "allelic_bams"
        conda: CONDA_SHARED_ENV
        shell:
            "SNPsplit {params.pairedEnd}"
            " -o {params.outdir} --snp_file {input.snp} {input.bam} 2> {log}"

elif aligner == "STAR" or aligner == "EXTERNAL_BAM":
    rule snp_split:
        input:
            snp = SNPFile,
            bam = aligner+"/{sample}.bam"
        output:
            targetbam = expand("allelic_bams/{{sample}}.{suffix}.bam", suffix = ['allele_flagged', 'genome1', 'genome2', 'unassigned']),
            tempbam = temp(aligner+"/{sample}.sortedByName.bam"),
            rep1 = "allelic_bams/{sample}.SNPsplit_report.yaml",
            rep2 = "allelic_bams/{sample}.SNPsplit_sort.yaml"
        log: "allelic_bams/logs/{sample}.snp_split.log"
        params:
            pairedEnd = '--paired' if pairedEnd else '',
            outdir = "allelic_bams"
        conda: CONDA_SHARED_ENV
        shell:
            "SNPsplit {params.pairedEnd}"
            " -o {params.outdir} --snp_file {input.snp} {input.bam} 2> {log}"

# move the allele-specific bams to another folder
#if aligner == "Bowtie2":
#    rule movebams:
#        input:
#            "allelic_bams/{sample}.filtered.{suffix}.bam"
#        params:
#            otherFile = "filtered_bam/{sample}.filtered.sortedByName.bam"
#        output:
#            temp("allelic_bams/{sample}.{suffix}.unsorted.bam")
#        shell:
#            "mv {input} {output} && \
#            if [ -e {params.otherFile} ]; then rm {params.otherFile}; fi"
#else:
#    rule movebams:
#        input:
#            "allelic_bams/{sample}.{suffix}.bam"
#        params:
#            otherFile = aligner+"/{sample}.sortedByName.bam"
#        output:
#            temp("allelic_bams/{sample}.{suffix}.unsorted.bam")
#        shell:
#            "mv {input} {output} && \
#            if [ -e {params.otherFile} ]; then rm {params.otherFile}; fi"

# sort them
rule BAMsort_allelic:
    input: "allelic_bams/{sample}.filtered.{suffix}.bam" if aligner == "Bowtie2" else "allelic_bams/{sample}.{suffix}.bam"
    output:
        "allelic_bams/{sample}.{suffix}.sorted.bam"
    log: "allelic_bams/logs/{sample}.{suffix}.sort.log"
    threads:
        12
    params:
        tempDir = tempDir
    conda: CONDA_SHARED_ENV
    shell: """
        TMPDIR={params.tempDir}
        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
        samtools sort -@ {threads} -T $MYTEMP -O bam -o {output} {input} 2> {log};
        rm -rf $MYTEMP
        """

# index the sorted files
rule BAMindex_allelic:
    input:
        "allelic_bams/{sample}.{suffix}.sorted.bam"
    output:
        "allelic_bams/{sample}.{suffix}.sorted.bam.bai"
    log: "allelic_bams/logs/{sample}.{suffix}.index.log"
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input} 2> {log}"
