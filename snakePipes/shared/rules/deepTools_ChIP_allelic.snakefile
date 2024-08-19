### deepTools bamCompare log2ratio #######################################################


rule bamCompare_log2_genome1:
    input:
        chip_bam = "allelic_bams/{chip_sample}.genome1.sorted.bam",
        chip_bai = "allelic_bams/{chip_sample}.genome1.sorted.bam.bai",
        control_bam = lambda wildcards: "allelic_bams/"+get_control(wildcards.chip_sample)+".genome1.sorted.bam",
        control_bai = lambda wildcards: "allelic_bams/"+get_control(wildcards.chip_sample)+".genome1.sorted.bam.bai"
    output:
        "deepTools_ChIP/bamCompare/allele_specific/{chip_sample}.genome1.log2ratio.over_{control_name}.bw"#, chip_sample = chip_samples, suffix = ['genome1', 'genome2'])
    conda: CONDA_SHARED_ENV
    params:
        bwBinSize = bwBinSize,
        ignoreForNorm = "--ignoreForNormalization " + ignoreForNormalization if ignoreForNormalization else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads " + str(fragmentLength),
        blacklist = "--blackListFileName " + blacklist_bed if blacklist_bed
                    else "",
        scaleFactors = " --scaleFactorsMethod readCount "
    benchmark:
        "deepTools_ChIP/.benchmark/bamCompare.log2ratio.{chip_sample}.{control_name}.genome1.benchmark"
    threads: lambda wildcards: 16 if 16<max_thread else max_thread
    shell: bamcompare_log2_cmd

rule bamCompare_log2_genome2:
    input:
        chip_bam = "allelic_bams/{chip_sample}.genome2.sorted.bam",
        chip_bai = "allelic_bams/{chip_sample}.genome2.sorted.bam.bai",
        control_bam = lambda wildcards: "allelic_bams/"+get_control(wildcards.chip_sample)+".genome2.sorted.bam",
        control_bai = lambda wildcards: "allelic_bams/"+get_control(wildcards.chip_sample)+".genome2.sorted.bam.bai"
    output:
        "deepTools_ChIP/bamCompare/allele_specific/{chip_sample}.genome2.log2ratio.over_{control_name}.bw"
    conda: CONDA_SHARED_ENV
    params:
        bwBinSize = bwBinSize,
        ignoreForNorm = "--ignoreForNormalization " + ignoreForNormalization if ignoreForNormalization else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads " + str(fragmentLength),
        blacklist = "--blackListFileName " + blacklist_bed if blacklist_bed
                    else "",
        scaleFactors = " --scaleFactorsMethod readCount "
    benchmark:
        "deepTools_ChIP/.benchmark/bamCompare.log2ratio.{chip_sample}.{control_name}.genome2.benchmark"
    threads: lambda wildcards: 16 if 16<max_thread else max_thread
    shell: bamcompare_log2_cmd

### deepTools plotEnrichment ###################################################

rule plotEnrichment_allelic:
    input:
        bams = expand("allelic_bams/{sample}.{suffix}.sorted.bam", sample = all_samples, suffix = ['genome1', 'genome2']),
        bais = expand("allelic_bams/{sample}.{suffix}.sorted.bam.bai", sample = all_samples, suffix = ['genome1', 'genome2'])
    output:
        png = "deepTools_ChIP/plotEnrichment/plotEnrichment.gene_features_allelic.png",
        tsv = "deepTools_ChIP/plotEnrichment/plotEnrichment.gene_features_allelic.tsv",
    conda: CONDA_SHARED_ENV
    params:
        genes_gtf = genes_gtf,
        labels = " --labels " + " ".join(expand("{sample}_{suffix}", sample = all_samples, suffix = ['genome1', 'genome2'])),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads " + str(fragmentLength)
    benchmark:
        "deepTools_ChIP/.benchmark/plotEnrichment_allelic.benchmark"
    threads: lambda wildcards: 24 if 24<max_thread else max_thread
    shell: plotEnrich_chip_cmd


### deepTools plotFingerprint (all files) ######################################

rule plotFingerprint_allelic:
    input:
        bams = expand("allelic_bams/{sample}.{suffix}.sorted.bam", sample = all_samples, suffix = ['genome1', 'genome2']),
        bais = expand("allelic_bams/{sample}.{suffix}.sorted.bam.bai", sample = all_samples, suffix = ['genome1', 'genome2'])
    output:
        metrics = "deepTools_ChIP/plotFingerprint/plotFingerprint.metrics_allelic.txt"
    conda: CONDA_SHARED_ENV
    params:
        labels = " --labels " + " ".join(expand("{sample}_{suffix}", sample = all_samples, suffix = ['genome1', 'genome2'])),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads " + str(fragmentLength),
        png = "--plotFile deepTools_ChIP/plotFingerprint/plotFingerprint_allelic.png" if (len(all_samples)<=20)
              else "",
        jsd = ""
    benchmark:
        "deepTools_ChIP/.benchmark/plotFingerprint_allelic.benchmark"
    threads: lambda wildcards: 24 if 24<max_thread else max_thread
    shell: plotFingerprint_cmd
