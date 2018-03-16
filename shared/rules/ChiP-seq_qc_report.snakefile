#####QC report for each sample########
rule qc_report_sample:
     input:
          duplication = "Sambamba/{sample}.dup.converted.tsv",
          peak = "MACS2/{sample}.filtered.BAM_peaks.qc.txt"
     output:
          tsv = "QC_report/{sample}.QCreport.tsv"
     run:
         shell("tail -n 1 {input.peak} > MACS2/{wildcards.sample}.temp"),
         shell("paste {input.duplication} MACS2/{wildcards.sample}.temp | column -s $'\t' -t | tee {output.tsv}"),
         shell("rm MACS2/{wildcards.sample}.temp")
########QC report for all the samples#########
rule qc_report_all:
     input:
         expand("QC_report/{sample}.QCreport.tsv",sample=all_samples)
     output:
         "QC_report/QC_report_all.tsv"
     run:
        header = (  "Sample\t"
                    "Total QC-passed reads(P)\t"
                    "Total QC-failed reads(F)\t"
                    "Duplicated P\t"
                    "Duplicated F\t"
                    "Mapped P\t"
                    "Mapped F\t"
                    "peak_count\t"
                    "FRiP\t"
                    "peak_genome_coverage"
                )
        shell('cat {input} | cat <(echo "'+header+'") - | tee {output}') 
