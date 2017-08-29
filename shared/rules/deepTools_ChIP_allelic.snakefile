### deepTools bamCompare log2ratio #######################################################

rule bamCompare_log2_genome1:
    input:
        chip_bam = "allelic_bams/{chip_sample}.genome1.sorted.bam",
        chip_bai = "allelic_bams/{chip_sample}.genome1.sorted.bam.bai",
        control_bam = lambda wildcards: "allelic_bams/"+get_control(wildcards.chip_sample)+".genome1.sorted.bam",
        control_bai = lambda wildcards: "allelic_bams/"+get_control(wildcards.chip_sample)+".genome1.sorted.bam.bai"
    output:
        "deepTools_ChIP/bamCompare/allele_specific/{chip_sample}.genome1.log2ratio.over_input.bw"#, chip_sample = chip_samples, suffix = ['genome1', 'genome2'])
    params:
        bw_binsize = bw_binsize,
        ignoreForNorm = "--ignoreForNormalization " + ignore_forNorm if ignore_forNorm else "",
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
    log:
        "deepTools_ChIP/logs/bamCompare.log2ratio.{chip_sample}.genome1.log"
    benchmark:
        "deepTools_ChIP/.benchmark/bamCompare.log2ratio.{chip_sample}.genome1.benchmark"
    threads: 16
    run:
        shell(bamcompare_log2_cmd())

rule bamCompare_log2_genome2:
    input:
        chip_bam = "allelic_bams/{chip_sample}.genome2.sorted.bam",
        chip_bai = "allelic_bams/{chip_sample}.genome2.sorted.bam.bai",
        control_bam = lambda wildcards: "allelic_bams/"+get_control(wildcards.chip_sample)+".genome2.sorted.bam",
        control_bai = lambda wildcards: "allelic_bams/"+get_control(wildcards.chip_sample)+".genome2.sorted.bam.bai"
    output:
        "deepTools_ChIP/bamCompare/allele_specific/{chip_sample}.genome2.log2ratio.over_input.bw"
    params:
        bw_binsize = bw_binsize,
        ignoreForNorm = "--ignoreForNormalization " + ignore_forNorm if ignore_forNorm else "",
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
    log:
        "deepTools_ChIP/logs/bamCompare.log2ratio.{chip_sample}.genome2.log"
    benchmark:
        "deepTools_ChIP/.benchmark/bamCompare.log2ratio.{chip_sample}.genome2.benchmark"
    threads: 16
    run:
        shell(bamcompare_log2_cmd())

### deepTools plotEnrichment ###################################################

rule plotEnrichment_allelic:
    input:
        bams = expand("allelic_bams/{sample}.{suffix}.sorted.bam", sample = all_samples, suffix = ['genome1', 'genome2']),
        bais = expand("allelic_bams/{sample}.{suffix}.sorted.bam.bai", sample = all_samples, suffix = ['genome1', 'genome2'])
    output:
        png = "deepTools_ChIP/plotEnrichment/plotEnrichment.gene_features_allelic.png",
        tsv = "deepTools_ChIP/plotEnrichment/plotEnrichment.gene_features_allelic.tsv",
    params:
        genes_gtf = genes_gtf,
        labels = " ".join(expand("{sample}_{suffix}", sample = all_samples, suffix = ['genome1', 'genome2'])),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length)
    log:
        "deepTools_ChIP/logs/plotEnrichment_allelic.log"
    benchmark:
        "deepTools_ChIP/.benchmark/plotEnrichment_allelic.benchmark"
    threads: 24
    run:
        shell(plotEnrich_chip_cmd())


### deepTools plotFingerprint (all files) ######################################

rule plotFingerprint_allelic:
    input:
        bams = expand("allelic_bams/{sample}.{suffix}.sorted.bam", sample = all_samples, suffix = ['genome1', 'genome2']),
        bais = expand("allelic_bams/{sample}.{suffix}.sorted.bam.bai", sample = all_samples, suffix = ['genome1', 'genome2'])
    output:
        metrics = "deepTools_ChIP/plotFingerprint/plotFingerprint.metrics_allelic.txt"
    params:
        labels = " ".join(all_samples),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length),
        png = "--plotFile deepTools_ChIP/plotFingerprint/plotFingerprint_allelic.png" if (len(all_samples)<=20)
              else "",
        jsd = ""
    log:
        "deepTools_ChIP/logs/plotFingerprint_allelic.log"
    benchmark:
        "deepTools_ChIP/.benchmark/plotFingerprint_allelic.benchmark"
    threads: 24
    run:
        shell(plotFingerprint_cmd())
