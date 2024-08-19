rule plotFingerprint:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample = samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample = samples)
    output:
        metrics = os.path.join(deeptools_ATAC, "plotFingerprint/plotFingerprint.metrics.txt")
    params:
        labels = " --labels " + " ".join(samples),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
        read_extension = "--extendReads",
        png = "--plotFile " + os.path.join(deeptools_ATAC, "plotFingerprint", "plotFingerprint.png") if (len(samples)<=20)
            else "",
        jsd = "--JSDsample filtered_bam/{}.filtered.bam".format(samples[0]) if (len(samples)>0)
            else ""
    benchmark:
        os.path.join(deeptools_ATAC, ".benchmark/plotFingerprint.benchmark")
    threads: lambda wildcards: 24 if 24<max_thread else max_thread
    conda: CONDA_SHARED_ENV
    shell: plotFingerprint_cmd


rule plotFingerprint_allelic:
    input:
        bams = expand("allelic_bams/{sample}.{suffix}.sorted.bam", sample = samples, suffix = ['genome1', 'genome2']),
        bais = expand("allelic_bams/{sample}.{suffix}.sorted.bam.bai", sample = samples, suffix = ['genome1', 'genome2'])
    output:
        metrics = os.path.join(deeptools_ATAC, "plotFingerprint", "plotFingerprint.metrics_allelic.txt")
    params:
        labels = " --labels " + " ".join(expand("{sample}_{suffix}", sample = samples, suffix = ['genome1', 'genome2'])),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed
                    else "",
        read_extension = "--extendReads",
        png = "--plotFile {}".format(os.path.join(deeptools_ATAC, "plotFingerprint", "plotFingerprint_allelic.png")) if (len(samples)<=20)
              else "",
        jsd = ""
    benchmark:
        os.path.join(deeptools_ATAC, ".benchmark/plotFingerprint_allelic.benchmark")
    threads: lambda wildcards: 24 if 24<max_thread else max_thread
    conda: CONDA_SHARED_ENV
    shell: plotFingerprint_cmd


# samtools, gawk
# bc from conda segfaults
rule MACS2_peak_qc:
    input:
        bam = os.path.join(short_bams, "{sample}.short.cleaned.bam"),
        bai = os.path.join(short_bams, "{sample}.short.cleaned.bam.bai"),
        xls = os.path.join(outdir_MACS2, '{sample}.filtered.short.BAM_peaks.xls')
    output:
        qc = os.path.join(outdir_ATACqc, "{sample}.filtered.BAM_peaks.qc.txt")
    params:
        peaks = os.path.join(outdir_MACS2, '{sample}.filtered.short.BAM_peaks.narrowPeak'),
        genome_index = genome_index
    benchmark:
        os.path.join(outdir_ATACqc, ".benchmark/ATAC_qc.{sample}.filtered.benchmark")
    conda: CONDA_SHARED_ENV
    shell: """
        # get the number of peaks
        peak_count=`cat {params.peaks} | wc -l`

        # get the number of mapped reads
        mapped_reads=`samtools view -c -F 4 {input.bam}`

        # calculate the number of alignments overlapping the peaks
        # exclude reads flagged as unmapped (unmapped reads will be reported when using -L)
        reads_in_peaks=`samtools view -c -F 4 -L {params.peaks} {input.bam}`

        # calculate Fraction of Reads In Peaks
        frip=`bc -l <<< "$reads_in_peaks/$mapped_reads"`

        # compute peak genome coverage
        peak_len=`awk '{{total+=$3-$2}}END{{print total}}' {params.peaks}`
        genome_size=`awk '{{total+=$3-$2}}END{{print total}}' {params.genome_index}`
        genomecov=`bc -l <<< "$peak_len/$genome_size"`

        # write peak-based QC metrics to output file
        printf "peak_count\tFRiP\tpeak_genome_coverage\n%d\t%5.3f\t%6.4f\n" $peak_count $frip $genomecov > {output.qc}
        """
