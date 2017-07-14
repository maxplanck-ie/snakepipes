if paired:
    rule featureCounts_allele:
        input:
            saf = "Annotation/genes.filtered.saf",
            bam = "allelic_bams/{sample}.allele_flagged.sorted.bam",
            allele1 = "allelic_bams/{sample}.genome1.sorted.bam",
            allele2 = "allelic_bams/{sample}.genome2.sorted.bam"
        output:
            'featureCounts/{sample}_allelic_counts.txt'
        params:
            libtype = library_type,
            opts = config["featurecounts_options"],
        log:
            "featureCounts/{sample}.log"
        threads: 8
        shell:
            feature_counts_path+"featureCounts"
            " -p -B {params.opts}"
            " -T {threads}"
            " -s {params.libtype}"
            " -F SAF -a {input.saf}"
            " -o {output}"
            " --tmpDir ${{TMPDIR}}"
            " {input.bam} {input.allele1} {input.allele2} &>> {log}"
else:
    rule featureCounts_allele:
        input:
            saf = "Annotation/genes.filtered.saf",
            bam = "allelic_bams/{sample}.allele_flagged.sorted.bam",
            allele1 = "allelic_bams/{sample}.genome1.sorted.bam",
            allele2 = "allelic_bams/{sample}.genome2.sorted.bam"
        output:
            'featureCounts/{sample}_allelic_counts.txt'
        params:
            libtype = library_type,
            opts = config["featurecounts_options"],
        log:
            "featureCounts/{sample}.log"
        threads: 8
        shell:
            feature_counts_path+"featureCounts "
            " {params.opts}"
            " -T {threads}"
            " -s {params.libtype}"
            " -F SAF -a {input.saf}"
            " -o {output}"
            " --tmpDir ${{TMPDIR}}"
            " {input.bam} {input.allele1} {input.allele2} &>> {log}"

rule merge_featureCounts:
    input:
        expand("featureCounts/{sample}_allelic_counts.txt", sample=samples)
    output:
        "featureCounts/allelic_counts.txt"
    shell:
        R_path + "Rscript "+os.path.join(maindir, "shared", "tools", "merge_featureCounts.R")+" {output} {input}"
