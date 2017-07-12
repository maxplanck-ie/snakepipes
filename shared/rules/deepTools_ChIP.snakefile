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
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
    log:
        "deepTools_ChIP/logs/bamCompare.log2ratio.{chip_sample}.filtered.log"
    benchmark:
        "deepTools_ChIP/.benchmark/bamCompare.log2ratio.{chip_sample}.filtered.benchmark"
    threads: 16
    shell:
        deepTools_path+"bamCompare "
        "-b1 {input.chip_bam} "
        "-b2 {input.control_bam} "
        "-o {output} "
        "--ratio log2 "
        "--scaleFactorsMethod readCount "
        "--binSize {params.bw_binsize} "
        "-p {threads} "
        "{params.read_extension} "
        "{params.blacklist} "
        "&> {log}"


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
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                         else "",
    log:
        "deepTools_ChIP/logs/bamCompare.subtract.{chip_sample}.filtered.log"
    benchmark:
        "deepTools_ChIP/.benchmark/bamCompare.subtract.{chip_sample}.filtered.benchmark"
    threads: 16
    shell:
        deepTools_path+"bamCompare "
        "-b1 {input.chip_bam} "
        "-b2 {input.control_bam} "
        "-o {output} "
        "--ratio subtract "
        "--scaleFactorsMethod readCount "
        "--normalizeTo1x {params.genome_size} "
        "--binSize {params.bw_binsize} "
        "-p {threads} "
        "{params.read_extension} "
        "{params.blacklist} "
        "&> {log}"


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
    shell:
        deepTools_path+"plotEnrichment "
        "-b {input.bams} "
        "--BED {params.genes_gtf} "
        "--plotFile {output.png} "
        "--labels {params.labels} "
        "--plotTitle 'Sigal enrichment (fraction of reads) without duplicates' "
        "--outRawCounts {output.tsv} "
        "--variableScales "
        "{params.blacklist} "
        "-p {threads} "
        "{params.read_extension} "
        "--ignoreDuplicates "
        "&> {log}"


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
        "{params.jsd} "
        "&> {log}"
