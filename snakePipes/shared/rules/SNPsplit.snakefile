
## get input bam depending on the mapping prog (use filtered bam in case of chip-seq data)
if config['aligner'] == "Bowtie2":
    rule snp_split:
        input:
            snp = config['SNPFile'],
            bam = "filtered_bam/{sample}.filtered.bam"
        output:
            targetbam = expand("allelic_bams/{{sample}}.filtered.{suffix}.bam", suffix = ['allele_flagged', 'genome1', 'genome2', 'unassigned']),
            tempbam = temp("filtered_bam/{sample}.filtered.sortedByName.bam")
        log: "allelic_bams/logs/{sample}.snp_split.log"
        params:
            pairedEnd = '--paired' if config['pairedEnd'] else '',
            outdir = "allelic_bams"
        conda: config['CONDA_SHARED_ENV']
        shell:
            "SNPsplit {params.pairedEnd}"
            " -o {params.outdir} --snp_file {input.snp} {input.bam} 2> {log}"

elif config['aligner'] == "STAR":
    rule snp_split:
        input:
            snp = config['SNPFile'],
            bam = config['aligner']+"/{sample}.bam"
        output:
            targetbam = expand("allelic_bams/{{sample}}.{suffix}.bam", suffix = ['allele_flagged', 'genome1', 'genome2', 'unassigned']),
            tempbam = temp(config['aligner']+"/{sample}.sortedByName.bam")
        log: "allelic_bams/logs/{sample}.snp_split.log"
        params:
            pairedEnd = '--paired' if config['pairedEnd'] else '',
            outdir = "allelic_bams"
        conda: config['CONDA_SHARED_ENV']
        shell:
            "SNPsplit {params.pairedEnd}"
            " -o {params.outdir} --snp_file {input.snp} {input.bam} 2> {log}"

# move the allele-specific bams to another folder
#if config['aligner'] == "Bowtie2":
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
#            otherFile = config['aligner']+"/{sample}.sortedByName.bam"
#        output:
#            temp("allelic_bams/{sample}.{suffix}.unsorted.bam")
#        shell:
#            "mv {input} {output} && \
#            if [ -e {params.otherFile} ]; then rm {params.otherFile}; fi"

# sort them
rule BAMsort_allelic:
    input: "allelic_bams/{sample}.filtered.{suffix}.bam" if config['aligner'] == "Bowtie2" else "allelic_bams/{sample}.{suffix}.bam"
    output:
        "allelic_bams/{sample}.{suffix}.sorted.bam"
    log: "allelic_bams/logs/{sample}.{suffix}.sort.log"
    threads:
        12
    params:
        tempDir = config['tempDir']
    conda: config['CONDA_SHARED_ENV']
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
    conda: config['CONDA_SHARED_ENV']
    shell: "samtools index {input} 2> {log}"
