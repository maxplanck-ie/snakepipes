deeptools_ATAC='deepTools_ATAC'
#
### TODO: Plot coverage from cutsites,based on the shifted the 5' end of
###         each fragment.
#
# rule bamCoverageCutsite:
#     input:
#         bam = "filtered_bam/{sample}.filtered.cutsites.bam",
#         bai = "filtered_bam/{sample}.filtered.cutsites.bam.bai"
#     output:
#         os.path.join(deeptools_ATAC, "{sample}.cutsites.bw")
#     # params:
#     #     pass
#     log:
#         pass
#     benchmark:
#         pass
#     shell:
#         pass
        # 1) convert bedpe to bam
        # 2) produce bw tracks
        # 3) remove bam

rule bamCoverageFragments:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        bai = "filtered_bam/{sample}.filtered.bam"
    output:
        os.path.join(deeptools_ATAC, "{sample}.fragments.bw" )
    params:
        bw_binsize = bw_binsize,
        genome_size = genome_size,
        read_extension = "--extendReads"
    threads: 16
    shell:
        deepTools_path+"bamCoverage "
        "-b {input.bam} "
        "-o {output} "
        "--binSize {params.bw_binsize} "
        "-p {threads} "
        "--normalizeTo1x {params.genome_size} "
        "{params.read_extension} "
        "&> {log}"



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
              else ""
        jsdsample="--JSDsample " + samples[0]
    log:
        os.path.join(deeptools_ATAC, "logs","plotFingerprint.log")
    benchmark:
        os.path.join(deeptools_ATAC, ".benchmark","plotFingerprint.benchmark")
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
        "{params.jsdsample} "
        "&> {log}"
