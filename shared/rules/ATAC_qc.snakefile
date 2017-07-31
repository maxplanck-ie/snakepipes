def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

MACS2_qc_folder = 'MACS2_qc/'
deeptools_ATAC = "deepTools_ATAC"

rule MACS2_peak_qc:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        aln_metrics = "Picard_qc/AlignmentSummaryMetrics/{sample}.alignment_summary_metrics.txt",
        xls = 'peaks_openChromatin/openchromatin_{sample}_peaks.xls'
    output:
        qc = MACS2_qc_folder + "{sample}.filtered.BAM_peaks.qc.txt"
    params:
        peaks = 'peaks_openChromatin/openchromatin_{sample}_peaks.narrowPeak',
        genome_index = genome_index
    log:
        MACS2_qc_folder + "logs/MACS2_peak_qc.{sample}.filtered.log"
    benchmark:
        MACS2_qc_folder + ".benchmark/MACS2_peak_qc.{sample}.filtered.benchmark"
    run:
        # get the number of peaks
        cmd = "cat "+params.peaks+" | wc -l"
        peak_count = int(subprocess.check_output( cmd, shell=True).decode())

        # get the number of mapped reads from Picard CollectAlignmentSummaryMetrics output
        cmd = "egrep '^PAIR|UNPAIRED' "+input.aln_metrics+" | cut -f 6"
        mapped_reads = int(subprocess.check_output( cmd, shell=True).decode())

        # calculate the number of alignments overlapping the peaks
        # exclude reads flagged as unmapped (unmapped reads will be reported when using -L)
        cmd = samtools_path+"samtools view -c -F 4 -L "+params.peaks+" "+input.bam
        reads_in_peaks = int(subprocess.check_output( cmd, shell=True).decode())

        # calculate Fraction of Reads In Peaks
        frip = reads_in_peaks / mapped_reads

        # compute peak genome coverage
        cmd = ("sort -k 1,1 "+params.peaks+" | "+
               bedtools_path+"genomeCoverageBed -i - -g "+params.genome_index+" | "+
               "grep -P '^genome\t1' | cut -f 5"
              )
        res=subprocess.check_output( cmd, shell=True).decode()
        genomecov=0
        if isFloat(res):
        	genomecov=float(res)

        # write peak-based QC metrics to output file
        with open(output.qc, "w") as f:
            f.write("peak_count\tFRiP\tpeak_genome_coverage\n"
                    "{:d}\t{:.3f}\t{:.4f}\n".format(
                    peak_count, frip, genomecov))

rule plotFingerprint:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample = samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample = samples)
    output:
        metrics = "deepTools_ATAC/plotFingerprint/plotFingerprint.metrics.txt"
    params:
        labels = " ".join(samples),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
        read_extension = "--extendReads",
        png = "--plotFile " + deeptools_ATAC + "/plotFingerprint/plotFingerprint.png" if (len(samples)<=20)
              else ""
    log:
        deeptools_ATAC + "/logs/plotFingerprint.log"
    benchmark:
        deeptools_ATAC + "/.benchmark/plotFingerprint.benchmark"
    threads: 24
    shell:
        deepTools_path+"plotFingerprint "
        "-b {input.bams} "
        "--labels {params.labels} "
        "--plotTitle 'Cumulative read counts per bin without duplicates' "
        "--ignoreDuplicates "
        "--outQualityMetrics {output.metrics} "
        "-p {threads} "
        "{params.blacklist} "
        "{params.png} "
        "{params.read_extension} "
        "&> {log}"
