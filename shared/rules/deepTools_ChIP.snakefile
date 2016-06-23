### deepTools bamCompare #######################################################

rule bamCompare:
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
                         else "--extendReads "+str(fragment_length)
    log:
        "deepTools_ChIP/logs/bamCompare.{chip_sample}.filtered.log"
    benchmark:
        "deepTools_ChIP/.benchmark/bamCompare.{chip_sample}.filtered.benchmark"
    threads: 16
    shell:
        deepTools_path+"bamCompare "
        "-b1 {input.chip_bam} "
        "-b2 {input.control_bam} "
        "-o {output} "
        "--ratio log2 "
        "--binSize {params.bw_binsize} "
        "-p {threads} "
        "{params.read_extension} "
        "&> {log}"


### deepTools plotEnrichment ###################################################

rule plotEnrichment:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample = all_samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample = all_samples)
    output:
        png = "deepTools_ChIP/plotEnrichment/signal_erichment.gene_features.png",
        tsv = "deepTools_ChIP/plotEnrichment/signal_erichment.gene_features.tsv"
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
# TODO: include blacklist parameter once the bug causing on error in plotEnrichment is fixed
#        "{params.blacklist} "
        "-p {threads} "
# TODO: include read extension parameter once the bug causing on error in plotEnrichment is fixed
#        "{params.read_extension} "
        "--ignoreDuplicates "
        "&> {log}"


### deepTools plotFingerprint ##################################################

rule plotFingerprint:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam", sample = all_samples),
        bais = expand("filtered_bam/{sample}.filtered.bam.bai", sample = all_samples)
    output:
        "deepTools_ChIP/plotFingerprint/fingerprint.png"
    params:
        labels = " ".join(all_samples),
        blacklist = "--blackListFileName "+blacklist_bed if blacklist_bed
                    else "",
        read_extension = "--extendReads" if paired
                         else "--extendReads "+str(fragment_length)
    log:
        "deepTools_ChIP/logs/plotFingerprint.log"
    benchmark:
        "deepTools_ChIP/.benchmark/plotFingerprint.benchmark"
    threads: 24
    shell:
        deepTools_path+"plotFingerprint "
        "-b {input.bams} "
        "--plotFile {output} "
        "--labels {params.labels} "
        "--plotTitle 'Cumulative read counts per bin without duplicates' "
        "{params.blacklist} "
        "-p {threads} "
        "{params.read_extension} "
        "--ignoreDuplicates "
        "&> {log}"
