
######### genomic contamination (IHEC metric)##############

rule GContamination_featureCounts:
    input:
        bams = aligner + "/{sample}.bam",
        gtf = str(extended_coding_regions_gtf or '')
    output:
        txt = temp("GenomicContamination/{sample}.featurecounts.txt"),
        summary = "GenomicContamination/{sample}.featurecounts.txt.summary"
    threads: 8
    conda: CONDA_RNASEQ_ENV
    shell:
        "featureCounts -T {threads} -a {input.gtf} -t transcript -o {output.txt} {input.bams}"

rule GContamination_featurecount_report:
    input:
        "GenomicContamination/{sample}.featurecounts.txt.summary"
    output:
        finaloutput = "GenomicContamination/{sample}.Gcontamination_report.tsv",
        temp = temp("GenomicContamination/{sample}.temp"),
        temp1 = temp("GenomicContamination/{sample}.temp1")
    shell:
        "SUM=$(cut -f2 {input} | tr '\n' '\t'| cut -f2,4,5 | awk '{{num = $1 + $2 + $3}} END {{print num}}');NUM=$(cut -f2 {input} | tr '\n' '\t'| cut -f5 | awk '{{num = $1}} END {{print num}}'); bc -l <<< $NUM/$SUM > {output.temp}; cut -f2 {input} | tr '\n' '\t'|cut -f2 | sed 's/^/{wildcards.sample}\t/' > {output.temp1}; paste -d'\t' {output.temp1} {output.temp} > {output.finaloutput}"

##
rule GContamination_featurecount_all_report:
    input:
        expand("GenomicContamination/{sample}.Gcontamination_report.tsv",sample=samples)
    output:
        report = "GenomicContamination/genomic_contamination_featurecount_report.tsv"
    shell: """
        echo "Sample\tMapped on transcripts\tGenomic contamination" > {output.report}
        sort -k1,1V {input} >> {output.report}
        """
