def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

rule MACS2_peak_qc:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        aln_metrics = "Picard_qc/AlignmentSummaryMetrics/{sample}.alignment_summary_metrics.txt",
        xls = os.path.join(outdir_MACS2, '{sample}_peaks.xls')
    output:
        qc = os.path.join(outdir_ATACqc, "{sample}.filtered.BAM_peaks.qc.txt")
    params:
        peaks = os.path.join(outdir_MACS2, '{sample}_peaks.narrowPeak'),
        genome_index = genome_index
    log:
        os.path.join(outdir_ATACqc, "logs/ATAC_qc.{sample}.filtered.log")
    benchmark:
        os.path.join(outdir_ATACqc, ".benchmark/ATAC_qc.{sample}.filtered.benchmark")
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
