### generate QC report for single sample ####################################

if paired:
    rule qc_report_sample:
        input:
            alignment_summary_metrics = "Picard_qc/AlignmentSummaryMetrics/{sample}.alignment_summary_metrics.txt",
            mark_duplicates_metrics = "Picard_qc/MarkDuplicates/{sample}.mark_duplicates_metrics.txt",
            insert_size_metrics = "Picard_qc/InsertSizeMetrics/{sample}.insert_size_metrics.txt",
            macs2_qc_txt = lambda wildcards: "MACS2/"+wildcards.sample+".filtered.BAM_peaks.qc.txt" if is_chip(wildcards.sample)
                else [],
            macs2_xls = lambda wildcards: "MACS2/"+wildcards.sample+".filtered.BAM_peaks.xls" if is_chip(wildcards.sample)
                else []
        output:
            tsv = temp("QC_report/{sample}.qc_report.tsv")
        log:
            "QC_report/logs/qc_report.{sample}.log"
        benchmark:
            "QC_report/.benchmark/qc_report.{sample}.benchmark"
        shell:
            "python " + os.path.join(workflow_tools, "sample_qc_report_PE.py") + " "
            "{input.alignment_summary_metrics} {input.mark_duplicates_metrics} {input.insert_size_metrics} {input.macs2_xls} {input.macs2_qc_txt} "
            ">{output} 2>{log} "
else:
    rule qc_report_sample:
        input:
            alignment_summary_metrics = "Picard_qc/AlignmentSummaryMetrics/{sample}.alignment_summary_metrics.txt",
            mark_duplicates_metrics = "Picard_qc/MarkDuplicates/{sample}.mark_duplicates_metrics.txt",
            macs2_qc_txt = lambda wildcards: "MACS2/"+wildcards.sample+".filtered.BAM_peaks.qc.txt" if is_chip(wildcards.sample)
                else [],
            macs2_xls = lambda wildcards: "MACS2/"+wildcards.sample+".filtered.BAM_peaks.xls" if is_chip(wildcards.sample)
                else []
        output:
            tsv = temp("QC_report/{sample}.qc_report.tsv")
        log:
            "QC_report/logs/qc_report.{sample}.log"
        benchmark:
            "QC_report/.benchmark/qc_report.{sample}.benchmark"
        shell:
            "python " + os.path.join(workflow_tools, "sample_qc_report_SE.py") + " "
            "{input.alignment_summary_metrics} {input.mark_duplicates_metrics} {input.macs2_xls} {input.macs2_qc_txt} "
            ">{output} 2>{log} "


### generate a single QC report for all samples ################################
rule qc_report_all:
    input:
        expand("QC_report/{sample}.qc_report.tsv", sample = all_samples)
    output:
        "QC_report/qc_report.all_samples.tsv"
    log:
        "QC_report/logs/qc_report.all_samples.log"
    benchmark:
        "QC_report/.benchmark/qc_report.all_samples.benchmark"
    run:
        if paired:
            header = (  "Sample\t"
                        "# Fragments\t"
                        "# Mapped fragments\t"
                        "Mapping rate\t"

                        "# Duplicated mapped fragments\t"
                        "Duplication rate\t"
                        "# Duplicate-free mapped fragments\t"
                        "Duplicate-free mapping rate\t"
                        "HQ read mapping rate (MAPQ>=20)\t"
                        "Singleton mapping rate\t"
                        "Median fragment size\t"

                        "# Peaks\t"
                        "FRiP\t"
                        "Peak genome coverage"
                    )
        else:
            header = (  "Sample\t"
                        "# Reads\t"
                        "# Mapped reads\t"
                        "Mapping rate\t"

                        "# Duplicated mapped reads\t"
                        "Duplication rate\t"
                        "# Duplicate-free mapped reads\t"
                        "Duplicate-free mapping rate\t"
                        "HQ read mapping rate (MAPQ>=20)\t"
                        "Fragment size (MACS2)\t"

                        "# Peaks\t"
                        "FRiP\t"
                        "Peak genome coverage"
                     )
        # output header and single sample QC sorted by sample name
        shell(
            'sort -k 1,1V {input} | '
            'cat <(echo "'+header+'") - '
            '> {output}'
        )
