#############convert sambamba report format########
rule convert_flagstat_output:
    input:
        "Sambamba/{sample}.markdup.txt"
    output:
        temp("Sambamba/{sample}.dup.converted.tsv")
    shell: """
        sed -n '1p;4p;5p' {input} | cut -d' ' -f1 | tr '\n' '\t' | sed 's/^/{wildcards.sample}\t/' | sed -e '$a\\' > {output}
        """

#######merge converted sambamba reports######
rule report_flagstat_all_data:
    input:
        expand("Sambamba/{sample}.dup.converted.tsv",sample=all_samples)
    output:
        "Sambamba/flagstat_report_all.tsv"
    shell: """
        echo -e 'sample\ttotal\tdup\tmapped' > {output}
        sort -k1,1V {input} >> {output}
        """

##########QC report for all the samples#########
rule qc_report_all:
    input:
        flagstat = "Sambamba/flagstat_report_all.tsv",
        metrics = "deepTools_ChIP/plotFingerprint/plotFingerprint.metrics.txt" if not useSpikeInForNorm else "split_deepTools_ChIP/plotFingerprint/plotFingerprint.metrics.txt"
    output:
        "QC_report/QC_report_all.tsv"
    shell: """
        awk 'NR == 1; NR > 1 {{print $0 | \"sort -k1,1V\"}}' {input.metrics} | cut -f4,8,10,12 | paste {input.flagstat} - > {output}
        """
