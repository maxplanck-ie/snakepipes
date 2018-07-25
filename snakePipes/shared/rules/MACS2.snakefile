import subprocess

# MACS2 should be called on already filtered, e.g. duplicate-free, BAM files
# for paired-end BAM files, Picard MarkDuplicates is fragment-based and
# therefore superior to MACS2 mate 1-based duplicate detection


### MACS2 peak calling #########################################################


if paired:
    rule MACS2:
        input:
            chip = "filtered_bam/{chip_sample}.filtered.bam",
            control =
                lambda wildcards: "filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam" if get_control(wildcards.chip_sample)
                else [],
            insert_size_metrics = "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv"
        output:
            peaks = "MACS2/{chip_sample}.filtered.BAM_peaks.xls",
            peaksPE = "MACS2/{chip_sample}.filtered.BAMPE_peaks.xls"
        params:
            fragment_length = lambda wildcards: cf.get_fragment_length("deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv", wildcards.chip_sample),
            genome_size = genome_size,
            # TODO: test BAMPE mode and activate BAMPE for paired-end data and BAM for single-end data
            # if results of BAMPE and BAM are in good agreement for paired-end data
            # does BAMPE mode really extends each read pair or does it only estimate a mean fragment size? the latter would be no advantage over BAM mode
            # format = "-f BAMPE" if paired else "-f BAM"
            broad_calling =
                lambda wildcards: "--broad" if is_broad(wildcards.chip_sample) else "",
            control_param =
                lambda wildcards: "-c filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam" if get_control(wildcards.chip_sample)
                else "",
        log:
            out = "MACS2/logs/MACS2.{chip_sample}.filtered.out",
            err = "MACS2/logs/MACS2.{chip_sample}.filtered.err"
        benchmark:
            "MACS2/.benchmark/MACS2.{chip_sample}.filtered.benchmark"
        conda: CONDA_CHIPSEQ_ENV
        shell: """
            macs2 callpeak -t {input.chip} \
                {params.control_param} -f BAM \
                -g {params.genome_size} --keep-dup all \
                --outdir MACS2 --name {wildcards.chip_sample}.filtered.BAM \
                --nomodel --extsize {params.fragment_length} {params.broad_calling} > {log.out} 2> {log.err}

            # also run MACS2 in paired-end mode BAMPE for comparison with single-end mode
            macs2 callpeak -t {input.chip} \
                {params.control_param} -f BAMPE \
                -g {params.genome_size} --keep-dup all \
                --outdir MACS2 --name {wildcards.chip_sample}.filtered.BAMPE \
                {params.broad_calling} > {log.out}.BAMPE 2> {log.err}.BAMPE
            """
else:
    rule MACS2:
        input:
            chip = "filtered_bam/{chip_sample}.filtered.bam",
            control =
                lambda wildcards: "filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam" if get_control(wildcards.chip_sample)
                else [],
        output:
            peaks = "MACS2/{chip_sample}.filtered.BAM_peaks.xls",
        params:
            fragment_length = fragment_length,
            genome_size = int(genome_size),
            broad_calling =
                lambda wildcards: "--broad" if is_broad(wildcards.chip_sample)
                else "",
            control_param =
                lambda wildcards: "-c filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam" if get_control(wildcards.chip_sample)
                else "",
        log:
            out = "MACS2/logs/MACS2.{chip_sample}.filtered.out",
            err = "MACS2/logs/MACS2.{chip_sample}.filtered.err"
        benchmark:
            "MACS2/.benchmark/MACS2.{chip_sample}.filtered.benchmark"
        conda: CONDA_CHIPSEQ_ENV
        shell: """
            macs2 callpeak -t {input.chip} {params.control_param} -f BAM -g {params.genome_size} --keep-dup all --outdir MACS2 \
                --name {wildcards.chip_sample}.filtered.BAM --nomodel --extsize {params.fragment_length} {params.broad_calling} > {log.out} 2> {log.err}
            """



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
            lambda wildcards: "MACS2/{}.filtered.BAM_peaks.broadPeak".format(wildcards.sample) if is_broad(wildcards.sample)
            else "MACS2/{}.filtered.BAM_peaks.narrowPeak".format(wildcards.sample),
        genome_index = genome_index
    benchmark:
        "MACS2/.benchmark/MACS2_peak_qc.{sample}.filtered.benchmark"
    conda: CONDA_CHIPSEQ_ENV
    shell: """
        # get the number of peaks
        peak_count=`wc -l < {params.peaks}`

        # get the number of mapped reads from Picard CollectAlignmentSummaryMetrics output
        mapped_reads=`egrep '^PAIR|UNPAIRED' {input.aln_metrics} | cut -f 6`

        # calculate the number of alignments overlapping the peaks
        # exclude reads flagged as unmapped (unmapped reads will be reported when using -L)
        reads_in_peaks=`samtools view -c -F 4 -L {params.peaks} {input.bam}`

        # calculate Fraction of Reads In Peaks
        frip=`bc -l <<< "$reads_in_peaks/$mapped_reads"`

        # compute peak genome coverage
        peak_len=`awk '{{total+=$3-$2}}END{{print total}}' {params.peaks}`
        genome_size=`awk '{{total+=$3-$2}}END{{print total}}' {params.peaks}`
        genomecov=`bc -l <<< "$peak_len/$genome_size"`

        # write peak-based QC metrics to output file
        printf "peak_count\tFRiP\tpeak_genome_coverage\n%d\t%5.3f\t%6.4f\n" $peak_count $frip $genomecov > {output.qc}
        """

# TODO
# add joined deepTools plotEnrichment call for all peaks and samples in one plot
