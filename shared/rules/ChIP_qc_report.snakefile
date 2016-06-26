### generate QC report for single sample ####################################

rule qc_report_sample:
    input:
        dup_metrics = "Picard_qc/MarkDuplicates/{sample}.mark_duplicates_metrics.txt",
        aln_metrics = "Picard_qc/AlignmentSummaryMetrics/{sample}.alignment_summary_metrics.txt",
        insert_size_metrics = "Picard_qc/InsertSizeMetrics/{sample}.insert_size_metrics.txt" if paired
            else [],
        peak_qc_file =
            lambda wildcards: "MACS2/"+wildcards.sample+".filtered.BAM_peaks.qc.txt" if is_chip(wildcards.sample)
            else [],
        peaks =
            lambda wildcards: "MACS2/"+wildcards.sample+".filtered.BAM_peaks.xls"  if is_chip(wildcards.sample)
            else []
    output:
        tsv = temp("QC_report/{sample}.qc_report.tsv")
    params:
        paired = paired,
        is_chip = lambda wildcards: is_chip(wildcards.sample)
    log:
        "QC_report/logs/qc_report.{sample}.log"
    benchmark:
        "QC_report/.benchmark/qc_report.{sample}.benchmark"
    run:
        # get alignment qc for read pairs/reads from Picard CollectAlignmentSummaryMetrics output
        cmd = "egrep '^PAIR|UNPAIRED' "+input.aln_metrics
        aln_qc = subprocess.check_output( cmd, shell=True).decode().split()
        # total reads (passing Illumina's filter)
        reads = int(aln_qc[2])
        # mapped reads
        mapped_reads = int(aln_qc[5])
        fmapped_reads = float(aln_qc[6])
        # fraction of high-quality mapped reads (MAPQ >=20)
        fhq_mapped_reads = int(aln_qc[8]) / reads

        if params.paired:
            read_pairs = int(reads / 2)
            # reads aligned in pairs
            mapped_reads_in_pairs = int(aln_qc[16])
            # aligned read pairs/fragments
            mapped_pairs = int(mapped_reads_in_pairs / 2)
            fmapped_pairs = mapped_pairs / read_pairs
            # singletons
            fmapped_singletons = (mapped_reads - mapped_reads_in_pairs) / reads

        # get duplicate qc from Picard MarkDuplicates output
        cmd = "egrep '^Unknown Library' "+input.dup_metrics
        dup_qc = subprocess.check_output( cmd, shell=True).decode().split("\t")
        if params.paired:
            dup_mapped_pairs = int(dup_qc[5])
            fdup_mapped_pairs = dup_mapped_pairs / mapped_pairs
            dupfree_mapped_pairs = mapped_pairs - dup_mapped_pairs
            # duplication-free mapping rate
            fdupfree_mapped_pairs = dupfree_mapped_pairs / read_pairs
        else:
            dup_mapped_reads = int(dup_qc[4])
            fdup_mapped_reads = dup_mapped_reads / mapped_reads
            dupfree_mapped_reads = mapped_reads - dup_mapped_reads
            # duplication-free mapping rate
            fdupfree_mapped_reads = dupfree_mapped_reads / reads

        # get fragment size from Picard CollectInsertSizeMetrics/MACS2 output
        if params.paired:
            fragment_size = get_fragment_length(input.insert_size_metrics)
        else:
            if params.is_chip:
                cmd = "grep '# d =' "+input.peaks
                fragment_size = int(subprocess.check_output( cmd, shell=True).decode().split()[3])
            else:
                fragment_size = "NA"

        # get peak qc from MACS2_peak_qc output
        if params.is_chip:
            with open(input.peak_qc_file, "r") as f:
                peak_qc = f.readlines()[1].strip()
        else:
            peak_qc = "NA\tNA\tNA"

        # write QC metrics to output file
        with open(output.tsv, "w") as f:
            if params.paired:
                f.write("{}\t{:10d}\t{:9d}\t{:.3f}\t{:9d}\t{:.3f}\t{:9d}\t{:.3f}\t{:.3f}\t{:.3f}\t{:9d}\t{}\n".format(
                            wildcards.sample,
                            read_pairs,
                            mapped_pairs,
                            fmapped_pairs,
                            dup_mapped_pairs,
                            fdup_mapped_pairs,
                            dupfree_mapped_pairs,
                            fdupfree_mapped_pairs,
                            fhq_mapped_reads,
                            fmapped_singletons,
                            fragment_size,
                            peak_qc )
                        )
            else:
                f.write("{}\t{:9d}\t{:9d}\t{:.3f}\t{:9d}\t{:.3f}\t{:9d}\t{:.3f}\t{:.3f}\t{:9d}\t{}\n".format(
                            wildcards.sample,
                            reads,
                            mapped_reads,
                            fmapped_reads,
                            dup_mapped_reads,
                            fdup_mapped_reads,
                            dupfree_mapped_reads,
                            fdupfree_mapped_reads,
                            fhq_mapped_reads,
                            fragment_size,
                            peak_qc )
                        )


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
