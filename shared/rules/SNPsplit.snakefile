# SNPsplit rule
if paired:
    SNPparam = '--paired'
else:
    SNPparam = ''

## get input bam depending on the mapping prog (use filtered bam in case of chip-seq data)
if mapping_prg == "Bowtie2":
    rule snp_split:
        input:
            snp = snp_file,
            bam = "filtered_bam/{sample}.filtered.bam"
        output:
            expand("allelic_bams/{{sample}}.filtered.{suffix}.bam", suffix = ['allele_flagged', 'genome1', 'genome2', 'unassigned'])
        params:
            paired = SNPparam,
            outdir = "allelic_bams"
        conda: CONDA_SHARED_ENV
        shell:
            "SNPsplit"
            " {params.paired} --samtools_path "+samtools_path+"samtools"
            " -o {params.outdir} --snp_file {input.snp} {input.bam}"
elif mapping_prg == "STAR":
    rule snp_split:
        input:
            snp = snp_file,
            bam = mapping_prg+"/{sample}.bam"
        output:
            expand("allelic_bams/{{sample}}.{suffix}.bam", suffix = ['allele_flagged', 'genome1', 'genome2', 'unassigned'])
        params:
            paired = SNPparam,
            outdir = "allelic_bams"
        conda: CONDA_SHARED_ENV
        shell:
            "SNPsplit"
            " {params.paired} --samtools_path "+samtools_path+"samtools"
            " -o {params.outdir} --snp_file {input.snp} {input.bam}"

# move the allele-specific bams to another folder
if mapping_prg == "Bowtie2":
    rule movebams:
        input:
            "allelic_bams/{sample}.filtered.{suffix}.bam"
        output:
            temp("allelic_bams/{sample}.{suffix}.unsorted.bam")
        shell:
            "mv {input} {output}"
else:
    rule movebams:
        input:
            "allelic_bams/{sample}.{suffix}.bam"
        output:
            temp("allelic_bams/{sample}.{suffix}.unsorted.bam")
        shell:
            "mv {input} {output}"

# sort them
rule BAMsort_allelic:
    input: "allelic_bams/{sample}.{suffix}.unsorted.bam"
    output:
        "allelic_bams/{sample}.{suffix}.sorted.bam"
    threads:
        12
    conda: CONDA_SHARED_ENV
    shell: "samtools sort -@ {threads} -T ${{TMPDIR}} -O bam -o {output} {input}"

# index the sorted files
rule BAMindex_allelic:
    input:
        "allelic_bams/{sample}.{suffix}.sorted.bam"
    output:
        "allelic_bams/{sample}.{suffix}.sorted.bam.bai"
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"
