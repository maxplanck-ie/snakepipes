### deepTools bamCompare subtract #######################################################

if bigWigType == "subtract" or bigWigType == "both":
    rule bamCompare_subtract:
        input:
            chip_bam = "split_bam/{chip_sample}_host.bam",
            chip_bai = "split_bam/{chip_sample}_host.bam.bai",
            control_bam = "split_bam/{control_name}_host.bam",
            control_bai = "split_bam/{control_name}_host.bam.bai"
        output:
            "split_deepTools_ChIP/bamCompare/{chip_sample}.host.subtract.{control_name}.bw"
        params:
            bwBinSize = bwBinSize,
            genome_size = genome_size,
            ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
            read_extension = "--extendReads" if pairedEnd else "--extendReads {}".format(fragmentLength),
            blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else ""
        log:
            out = "split_deepTools_ChIP/logs/bamCompare.subtract.{chip_sample}.host.subtract.{control_name}.out",
            err = "split_deepTools_ChIP/logs/bamCompare.subtract.{chip_sample}.host.subtract.{control_name}.err"
        benchmark:
            "split_deepTools_ChIP/.benchmark/bamCompare.subtract.{chip_sample}.host.subtract.{control_name}.benchmark"
        threads: 16
        conda: CONDA_SHARED_ENV
        shell: bamcompare_subtract_cmd

### deepTools bamCompare log2ratio #######################################################
if bigWigType == "log2ratio" or bigWigType == "both":
    rule bamCompare_log2:
        input:
            chip_bam = "split_bam/{chip_sample}_host.bam",
            chip_bai = "split_bam/{chip_sample}_host.bam.bai",
            control_bam = "split_bam/{control_name}_host.bam",
            control_bai = "split_bam/{control_name}_host.bam.bai",
        output:
            "split_deepTools_ChIP/bamCompare/{chip_sample}.host.log2ratio.over_{control_name}.bw"
        params:
            bwBinSize = bwBinSize,
            ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
            read_extension = "--extendReads" if pairedEnd else "--extendReads {}".format(fragmentLength),
            blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else ""
        log:
            out = "split_deepTools_ChIP/logs/bamCompare.log2ratio.{chip_sample}.{control_name}.host.out",
            err = "split_deepTools_ChIP/logs/bamCompare.log2ratio.{chip_sample}.{control_name}.host.err"
        benchmark:
            "split_deepTools_ChIP/.benchmark/bamCompare.log2ratio.{chip_sample}.{control_name}.host.benchmark"
        threads: 16
        conda: CONDA_SHARED_ENV
        shell: bamcompare_log2_cmd


### deepTools plotEnrichment ###################################################

rule plotEnrichment:
    input:
        bams = expand("split_bam/{sample}_host.bam", sample = all_samples),
        bais = expand("split_bam/{sample}_host.bam.bai", sample = all_samples)
    output:
        png = "split_deepTools_ChIP/plotEnrichment/plotEnrichment.gene_features.png",
        tsv = "split_deepTools_ChIP/plotEnrichment/plotEnrichment.gene_features.tsv",
    params:
        genes_gtf = genes_gtf,
        labels = " ".join(all_samples),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
        read_extension = "--extendReads" if pairedEnd else "--extendReads {}".format(fragmentLength)
    log:
        out = "split_deepTools_ChIP/logs/plotEnrichment.out",
        err = "split_deepTools_ChIP/logs/plotEnrichment.err"
    benchmark:
        "split_deepTools_ChIP/.benchmark/plotEnrichment.benchmark"
    threads: 24
    conda: CONDA_SHARED_ENV
    shell: plotEnrich_chip_cmd


### deepTools plotFingerprint (all files) ######################################

rule plotFingerprint:
    input:
        bams = expand("split_bam/{sample}_host.bam", sample = all_samples),
        bais = expand("split_bam/{sample}_host.bam.bai", sample = all_samples)
    output:
        metrics = "split_deepTools_ChIP/plotFingerprint/plotFingerprint.metrics.txt"
    params:
        labels = " --smartLabels ",
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
        read_extension = "--extendReads" if pairedEnd else "--extendReads {}".format(fragmentLength),
        png = "--plotFile deepTools_ChIP/plotFingerprint/plotFingerprint.png" if (len(samples)<=20)
              else "",
        jsd = "--JSDsample split_bam/{}.host.bam".format(control_samples[0]) if (len(control_samples)>0)
            else ""
    log:
        out = "split_deepTools_ChIP/logs/plotFingerprint.out",
        err = "split_deepTools_ChIP/logs/plotFingerprint.err"
    benchmark:
        "split_deepTools_ChIP/.benchmark/plotFingerprint.benchmark"
    threads: 24
    conda: CONDA_SHARED_ENV
    shell: plotFingerprint_cmd
