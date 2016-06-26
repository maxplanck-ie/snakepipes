import subprocess


### MACS2 peak calling #########################################################

rule MACS2:
    input:
        chip = "filtered_bam/{chip_sample}.filtered.bam",
        control =
            lambda wildcards: "filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam" if get_control(wildcards.chip_sample)
            else [],
        insert_size_metrics =
            "Picard_qc/InsertSizeMetrics/{chip_sample}.insert_size_metrics.txt" if paired
            else []
    output:
        peaks = "MACS2/{chip_sample}.filtered.BAM_peaks.xls",
        peaksPE = "MACS2/{chip_sample}.filtered.BAMPE_peaks.xls" if paired
                  else []
    params:
        fragment_length = fragment_length,
        paired = paired,
        genome_size = int(genome_size),
        # TODO: test BAMPE mode and activate BAMPE for paired-end data and BAMPE for single-end data
        # if results of BAMPE and BAM are in good agreement for paired-end data
        # does BAMPE mode really extends each read pair or does it only estimate a mean fragment size? the latter would be no advantage over BAM mode
        # format = "-f BAMPE" if paired else "-f BAM"
        broad_calling =
            lambda wildcards: "--broad" if is_broad(wildcards.chip_sample)
            else "",
        control_param =
            lambda wildcards: "-c filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam" if get_control(wildcards.chip_sample)
            else "",
    log:
        "MACS2/logs/MACS2.{chip_sample}.filtered.log"
    benchmark:
        "MACS2/.benchmark/MACS2.{chip_sample}.filtered.benchmark"
    run:
        if params.paired:
            median_fragment_length = get_fragment_length(input.insert_size_metrics)
        else:
            median_fragment_length = params.fragment_length
        model = "--nomodel --extsize "+str(median_fragment_length) if params.paired else "--extsize "+str(median_fragment_length)
        shell(
            macs2_path+"macs2 callpeak "
            "-t {input.chip} "
            "{params.control_param} "
            "-f BAM "
            "-g {params.genome_size} "
            # MACS2 should be called on already filtered, e.g. duplicate-free, BAM files
            # for paired-end BAM files, Picard MarkDuplicates is fragment-based and
            # therefore superior to MACS2 mate 1-based duplicate detection
            "--keep-dup all "
            "--outdir MACS2 "
            "--name {wildcards.chip_sample}.filtered.BAM "
            +model+" "
            "{params.broad_calling} "
            "&> {log}"
        )
        # also run MACS2 in paired-end mode BAMPE for comparison with single-end mode
        if params.paired:
            shell(
                macs2_path+"macs2 callpeak "
                "-t {input.chip} "
                "{params.control_param} "
                "-f BAMPE "
                "-g {params.genome_size} "
                # MACS2 should be called on already filtered, e.g. duplicate-free, BAM files
                "--keep-dup all "
                "--outdir MACS2 "
                "--name {wildcards.chip_sample}.filtered.BAMPE "
                "{params.broad_calling} "
                "&> {log}.BAMPE"
            )


### MACS2 peak quality control #################################################

rule MACS2_peak_qc:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        aln_metrics = "Picard_qc/AlignmentSummaryMetrics/{sample}.alignment_summary_metrics.txt",
        xls = "MACS2/{sample}.filtered.BAM_peaks.xls"
    output:
        qc = "MACS2/{sample}.filtered.BAM_peaks.qc.txt"
    params:
        peaks =
            lambda wildcards: "MACS2/{sample}.filtered.BAM_peaks.broadPeak" if is_broad(wildcards.sample)
            else "MACS2/{sample}.filtered.BAM_peaks.narrowPeak",
        genome_index = genome_index
    log:
        "MACS2/logs/MACS2_peak_qc.{sample}.filtered.log"
    benchmark:
        "MACS2/.benchmark/MACS2_peak_qc.{sample}.filtered.benchmark"
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
               "grep -P 'genome\t1' | cut -f 5"
              )
        genomecov = float(subprocess.check_output( cmd, shell=True).decode())

        # write peak-based QC metrics to output file
        with open(output.qc, "w") as f:
            f.write("peak_count\tFRiP\tpeak_genome_coverage\n"
                    "{:d}\t{:.3f}\t{:.4f}\n".format(
                    peak_count, frip, genomecov))

# TODO
# add joined deepTools plotEnrichment call for all peaks and samples in one plot
