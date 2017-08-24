### deepTools bamCompare subtract #######################################################

rule bamCompare_subtract:
    input:
        chip_bam = "filtered_bam/{chip_sample}.filtered.bam",
        chip_bai = "filtered_bam/{chip_sample}.filtered.bam.bai",
        control_bam = lambda wildcards: "filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam",
        control_bai = lambda wildcards: "filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam.bai",
    output:
        "deepTools_ChIP/bamCompare/{chip_sample}.filtered.subtract.input.bw"
    params:
        bw_binsize = bw_binsize,
        genome_size = genome_size,
        ignoreForNorm = "--ignoreForNormalization " + ignore_forNorm,
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                         else "",
    log:
        "deepTools_ChIP/logs/bamCompare.subtract.{chip_sample}.filtered.log"
    benchmark:
        "deepTools_ChIP/.benchmark/bamCompare.subtract.{chip_sample}.filtered.benchmark"
    threads: 16
    run:
        shell(bamcompare_subtract_cmd())

### deepTools bamCompare log2ratio #######################################################

rule bamCompare_log2:
    input:
        chip_bam = "filtered_bam/{chip_sample}.filtered.bam",
        chip_bai = "filtered_bam/{chip_sample}.filtered.bam.bai",
        control_bam = lambda wildcards: "filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam",
        control_bai = lambda wildcards: "filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam.bai",
    output:
        "deepTools_ChIP/bamCompare/{chip_sample}.filtered.log2ratio.over_input.bw"
    params:
        bw_binsize = bw_binsize,
        ignoreForNorm = "--ignoreForNormalization " + ignore_forNorm,
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
    log:
        "deepTools_ChIP/logs/bamCompare.log2ratio.{chip_sample}.filtered.log"
    benchmark:
        "deepTools_ChIP/.benchmark/bamCompare.log2ratio.{chip_sample}.filtered.benchmark"
    threads: 16
    run:
        shell(bamcompare_log2_cmd())


### deepTools plotEnrichment ###################################################

rule plotEnrichment:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample = all_samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample = all_samples)
    output:
        png = "deepTools_ChIP/plotEnrichment/plotEnrichment.gene_features.png",
        tsv = "deepTools_ChIP/plotEnrichment/plotEnrichment.gene_features.tsv",
    params:
        genes_gtf = genes_gtf,
        labels = " ".join(all_samples),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length)
    log:
        "deepTools_ChIP/logs/plotEnrichment.log"
    benchmark:
        "deepTools_ChIP/.benchmark/plotEnrichment.benchmark"
    threads: 24
    run:
        shell(plotEnrich_chip_cmd())


### deepTools plotFingerprint (all files) ######################################

rule plotFingerprint:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample = all_samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample = all_samples)
    output:
        metrics = "deepTools_ChIP/plotFingerprint/plotFingerprint.metrics.txt"
    params:
        labels = " ".join(all_samples),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length),
        png = "--plotFile deepTools_ChIP/plotFingerprint/plotFingerprint.png" if (len(all_samples)<=20)
              else "",
        jsd = "--JSDsample filtered_bam/"+control_samples[0]+".filtered.bam" if (len(control_samples)>0)
            else ""
    log:
        "deepTools_ChIP/logs/plotFingerprint.log"
    benchmark:
        "deepTools_ChIP/.benchmark/plotFingerprint.benchmark"
    threads: 24
    run:
        shell(plotFingerprint_cmd())
