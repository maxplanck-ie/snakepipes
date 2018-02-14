
######### genomic contamination (IHEC metric)##############

#rule gtf2bed:
#      input:
#          gtf = extended_coding_regions_gtf
#      output:
#          bed = "extendedRegion.bed"
#      shell:
#          "cat {input.gtf} | awk '{{OFS=\"\t\"; print $1,$4,$5,$2,$6,$7}}' > {output.bed}"
#
#rule GContamination:
#      input:
#          bam = mapping_prg+"/{sample}.bam",
#          bed ="extendedRegion.bed"
#      output:
#          "MappedtoIntergenicRegions/{sample}.intergenic.bam"
#      shell:
#          bedtools_path+"bedtools intersect -abam {input.bam} -b {input.bed} -v -s > {output}"
#

rule GContamination_featureCounts:
            input:
                bams = mapping_prg+"/{sample}.bam",
                gtf = extended_coding_regions_gtf
            output:
                txt = "GenomicContamination/{sample}.featurecounts.txt",
                summary = "GenomicContamination/{sample}.featurecounts.txt.summary"
            log:
                "GenomicContamination/featurecounts.log"
            threads: 8
            shell:
                feature_counts_path+"featureCounts -T {threads} -a {input.gtf} -t transcript -o {output.txt} {input.bams} &>> {log}"

rule GContamination_featurecount_report:
           input:
                "GenomicContamination/{sample}.featurecounts.txt.summary"
           output:
               "GenomicContamination/{sample}.Gcontamination_report.tsv"
           run:
#               "cols=2; for((i=1;i<=$cols;i++)); do  awk '{{print $'$i'}}' {input} | tr '\n' ' '; echo;  done > {output.temp}"
               shell("cat {input} | cut -f2 | tail -n +2 | tr '\n' '\t' > GenomicContamination/{wildcards.sample}.output.temp"),
               shell("cat GenomicContamination/{wildcards.sample}.output.temp | cut -f1,5 | sed 's/^/{wildcards.sample}\t/' > {output}"),
               shell("rm GenomicContamination/{wildcards.sample}.output.temp")
rule GContamination_featurecount_all_report:
         input:
             expand("GenomicContamination/{sample}.Gcontamination_report.tsv",sample=samples)
         output:
             report = "GenomicContamination/genomic_contamination_featurecount_report.tsv"
         run:
             header = (  "Sample\t"
                         "Mapped on transcripts\t"
                         "Genomics contamination"
                   )
             shell('cat {input} | cat <(echo "'+header+'") - | tee {output.report}')

