##########QC report for all the samples#########
rule qc_report_all:
     input:
         flagstat = "Sambamba/flagstat_report_all.tsv",
         metrics = "deepTools_ChIP/plotFingerprint/plotFingerprint.metrics.txt"
     output:
         "QC_report/QC_report_all.tsv"
     run:
        shell("sort -k 1 -1v {input.flagstat} -o {input.flagstat}") ,
        shell("sort -k 1 -1v {input.metrics} -o {input.metrics}"),
        shell("cat {input.metrics} | cut -k4,8,10,12 | paste {input.flagstat} {input.metrics} > {output})")

