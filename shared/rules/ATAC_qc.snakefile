

def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

rule plotFingerprint:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample = samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample = samples)
    output:
        metrics = os.path.join(deeptools_ATAC, "plotFingerprint/plotFingerprint.metrics.txt")
    params:
        labels = " ".join(samples),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
        read_extension = "--extendReads",
        png = "--plotFile " +os.path.join(deeptools_ATAC ,
         "plotFingerprint","plotFingerprint.png") if (len(samples)<=20)
            else "",
        jsd = "--JSDsample filtered_bam/"+samples[0]+".filtered.bam" if (len(samples)>0)
            else ""
    log:
        os.path.join(deeptools_ATAC, "logs","plotFingerprint.log")
    benchmark:
        os.path.join(deeptools_ATAC, ".benchmark","plotFingerprint.benchmark")
    threads: 24
    conda: CONDA_SHARED_ENV
    shell: plotFingerprint_cmd


rule plotFingerprint_allelic:
    input:
        bams = expand("allelic_bams/{sample}.{suffix}.sorted.bam", sample = samples, suffix = ['genome1', 'genome2']),
        bais = expand("allelic_bams/{sample}.{suffix}.sorted.bam.bai", sample = samples, suffix = ['genome1', 'genome2'])
    output:
        metrics = os.path.join(deeptools_ATAC, "plotFingerprint/plotFingerprint.metrics_allelic.txt")
    params:
        labels = " ".join(expand("{sample}_{suffix}", sample = samples, suffix = ['genome1', 'genome2'])),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
        read_extension = "--extendReads",
        png = "--plotFile " +os.path.join(deeptools_ATAC ,
         "plotFingerprint","plotFingerprint_allelic.png") if (len(samples)<=20)
            else "",
        jsd = ""
    log:
        os.path.join(deeptools_ATAC, "logs/plotFingerprint_allelic.log")
    benchmark:
        os.path.join(deeptools_ATAC, ".benchmark/plotFingerprint_allelic.benchmark")
    threads: 24
    conda: CONDA_SHARED_ENV
    shell: plotFingerprint_cmd


# samtools, bedtools
rule MACS2_peak_qc:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        aln_metrics = "Picard_qc/AlignmentSummaryMetrics/{sample}.alignment_summary_metrics.txt",
        xls = os.path.join(outdir_MACS2, '{sample}.filtered.BAM_peaks.xls')
    output:
        qc = os.path.join(outdir_ATACqc, "{sample}.filtered.BAM_peaks.qc.txt")
    params:
        peaks = os.path.join(outdir_MACS2, '{sample}.filtered.BAM_peaks.narrowPeak'),
        genome_index = genome_index
    log:
        os.path.join(outdir_ATACqc, "logs/ATAC_qc.{sample}.filtered.log")
    benchmark:
        os.path.join(outdir_ATACqc, ".benchmark/ATAC_qc.{sample}.filtered.benchmark")
    conda: CONDA_ATAC_ENV
    shell: """
        # get the number of peaks
        peak_count=`cat {params.peaks} | wc -l`

        # get the number of mapped reads from Picard CollectAlignmentSummaryMetrics output
        mapped_reads=`egrep '^PAIR|UNPAIRED' {input.aln_metrics} | cut -f 6`

        # calculate the number of alignments overlapping the peaks
        # exclude reads flagged as unmapped (unmapped reads will be reported when using -L)
        reads_in_peaks=`samtools view -c -F 4 -L {params.peaks} {input.bam}`

        # calculate Fraction of Reads In Peaks
        frip=`bc -l <<< "$reads_in_peaks/$mapped_reds"`

        # compute peak genome coverage
        genomecov=`sort -k 1,1 {params.peaks} | genomeCoverageBed -i - -g {params.genome_index} | grep -P '^genome\t1' | cut -f 5`
        #genomecov=0
        #if isFloat(res):
        #	genomecov=float(res)

        # write peak-based QC metrics to output file
        printf "peak_count\tFRiP\tpeak_genome_coverage\n%d\t%5.3f\t%6.4f\n" ${peak_count} ${frip} ${genomecov} > {output.qc}
        """
