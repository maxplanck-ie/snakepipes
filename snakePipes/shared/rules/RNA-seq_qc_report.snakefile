########convert sambamba report format########
rule convert_flagstat_output:
      input:
         "Sambamba/{sample}.markdup.txt"
      output:
         temp("Sambamba/{sample}.dup.converted.tsv")
      shell:
         "sed -n '1p;4p;5p' {input} | cut -d' ' -f1 | tr '\n' '\t' | sed 's/^/{wildcards.sample}\t/' | sed -e '$a\\' > {output}"

#######merge converted sambamba reports######
rule report_flagstat_all_data:
      input:
         expand("Sambamba/{sample}.dup.converted.tsv",sample=samples)
      output:
         "Sambamba/flagstat_report_all.tsv"
      shell:
         "sort -k1,1V {input} | cat <( echo -e 'sample\ttotal\tdup\tmapped') - > {output}"

##########QC report for all the samples#########
if dnaContam:
  rule qc_report_all:
        input:
            flagstat = "Sambamba/flagstat_report_all.tsv",
            IHECmetrics = "GenomicContamination/genomic_contamination_featurecount_report.tsv"
        output:
            "QC_report/QC_report_all.tsv"
        shell:
            "cut -f2,3 {input.IHECmetrics} | paste {input.flagstat} - > {output}"
else:
  rule qc_report_all:
        input:
            flagstat = "Sambamba/flagstat_report_all.tsv"
        output:
            "QC_report/QC_report_all.tsv"
        shell:
            "cp {input.flagstat} {output}"
