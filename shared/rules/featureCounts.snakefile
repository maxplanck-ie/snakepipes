if paired:
    rule featureCounts:
        input:
            bam = "HISAT2/{sample}.bam",
            saf = "Annotation/genes.filtered.saf",
        output:
            "featureCounts/{sample}/counts.txt"
        params:
            libtype = library_type,
            opts = config["featurecounts_options"],
        log:
            "featureCounts/{sample}/featureCounts.log"
        threads: 4
        shell:
            feature_counts_path+"featureCounts "
            "-p -B {params.opts} "
            "-T {threads} "
            "-s {params.libtype} "
            "-F SAF -a {input.saf} "
            "-o {output} "
            "{input.bam} &> {log} "
else:
    rule featureCounts:
        input:
            bam = "HISAT2/{sample}.bam",
            saf = "Annotation/genes.filtered.saf",
        output:
            "featureCounts/{sample}/counts.txt"
        params:
            libtype = library_type,
            opts = config["featurecounts_options"],
        log:
            "featureCounts/{sample}/featureCounts.log"
        threads: 4
        shell:
            feature_counts_path+"featureCounts "
            "{params.opts} "
            "-T {threads} "
            "-s {params.libtype} "
            "-F SAF -a {input.saf} "
            "-o {output} "
            "{input.bam} &> {log} "
