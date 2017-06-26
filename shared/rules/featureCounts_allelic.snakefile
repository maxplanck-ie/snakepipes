if paired:
    rule featureCounts_allele:
        input:
            saf = "Annotation/genes.filtered.saf",
            bam = mapping_prg+"/{sample}_nsorted.allele_flagged.bam",
            allele1 = mapping_prg+"/{sample}_nsorted.genome1.bam",
            allele2 = mapping_prg+"/{sample}_nsorted.genome2.bam"
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
            bam = mapping_prg+"/{sample}_nsorted.allele_flagged.bam",
            allele1 = mapping_prg+"/{sample}_nsorted.genome1.bam",
            allele2 = mapping_prg+"/{sample}_nsorted.genome2.bam"
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
        R_path + " Rscript "+os.path.join(maindir, "shared", "tools", "merge_featureCounts.R")+" {output} {input}"
