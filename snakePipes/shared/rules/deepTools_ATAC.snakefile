### deepTools bamCompare subtract #######################################################

rule bamCompare_subtract:
    input:
        chip_bam = os.path.join(outdir_MACS2, "{sample}.short.cleaned.bam"),
        chip_bai = os.path.join(outdir_MACS2, "{sample}.short.cleaned.bam.bai")
    output:
        "deepTools_ATAC/bamCompare/{sample}.filtered.bw"
    params:
        bw_binsize = bw_binsize,
        genome_size = genome_size,
        ignoreForNorm = "--ignoreForNormalization {}".format(ignore_forNorm) if ignore_forNorm else "",
        read_extension = "--extendReads" if paired else "--extendReads {}".format(fragment_length),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else ""
    log:
        out = "deepTools_ATAC/logs/bamCompare.{sample}.filtered.out",
        err = "deepTools_ATAC/logs/bamCompare.{sample}.filtered.out"
    benchmark:
        "deepTools_ATAC/.benchmark/deepTools_ATAC/logs/bamCompare.{sample}.filtered.benchmark"
    threads: 16
    conda: CONDA_SHARED_ENV
    shell: bamcov_cmd

