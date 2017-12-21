deeptools_ATAC='deepTools_ATAC'

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
    run:
        shell(plotFingerprint_cmd())
