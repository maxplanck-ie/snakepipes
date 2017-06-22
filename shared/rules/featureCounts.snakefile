if paired:
    rule featureCounts:
        input:
            bam = mapping_prg+"/{sample}.bam",
            saf = "Annotation/genes.filtered.saf",
        output:
            "featureCounts/{sample}.counts.txt"
        params:
            libtype = library_type,
            opts = config["featurecounts_options"],
        log:
            "featureCounts/{sample}.log"
        threads: 8
        shell:
            feature_counts_path+"featureCounts "
            "-p -B {params.opts} "
            "-T {threads} "
            "-s {params.libtype} "
            "-F SAF -a {input.saf} "
            "-o {output} "
            "--tmpDir ${{TMPDIR}} "
            "{input.bam} &>> {log} "
else:
    rule featureCounts:
        input:
            bam = mapping_prg+"/{sample}.bam",
            saf = "Annotation/genes.filtered.saf",
        output:
            "featureCounts/{sample}.counts.txt"
        params:
            libtype = library_type,
            opts = config["featurecounts_options"],
        log:
            "featureCounts/{sample}.log"
        threads: 8
        shell:
            feature_counts_path+"featureCounts "
            "{params.opts} "
            "-T {threads} "
            "-s {params.libtype} "
            "-F SAF -a {input.saf} "
            "-o {output} "
            "--tmpDir ${{TMPDIR}} "
            "{input.bam} &>> {log} "


rule merge_featureCounts:
    input:
        expand("featureCounts/{sample}.counts.txt", sample=samples)
    output:
        "featureCounts/counts.tsv"
    shell:
        R_path+"Rscript "+os.path.join(maindir, "shared", "tools", "merge_featureCounts.R")+" {output} {input}"
