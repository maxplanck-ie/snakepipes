### deepTools bamCompare subtract #######################################################

rule bamcoverage_short_cleaned:
    input:
        bam = os.path.join(short_bams, "{sample}.short.cleaned.bam"),
        bai = os.path.join(short_bams, "{sample}.short.cleaned.bam.bai")
    output:
        "deepTools_ATAC/bamCompare/{sample}.filtered.bw"
    params:
        bwBinSize = bwBinSize,
        genome_size = genome_size,
        ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
        read_extension = "--extendReads" if pairedEnd else "--extendReads {}".format(fragmentLength),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else ""
    benchmark:
        "deepTools_ATAC/.benchmark/deepTools_ATAC/logs/bamCompare.{sample}.filtered.benchmark"
    threads: lambda wildcards: 16 if 16<max_thread else max_thread
    conda: CONDA_SHARED_ENV
    shell: bamcov_cmd
