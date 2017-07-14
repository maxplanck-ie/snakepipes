
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
            expand("filtered_bam/{{sample}}.filtered.{suffix}.bam", suffix = ['allele_flagged', 'genome1', 'genome2', 'unassigned'])
        params:
            paired = SNPparam
        log:
            "allelic_bams/{sample}_SNPsplit.log"
        shell:
            SNPsplit_path + "SNPsplit"
            " {params.paired} --samtools_path "+samtools_path+"samtools"
            " --snp_file {input.snp} {input.bam} 2> {log}"
elif mapping_prg == "STAR":
    rule snp_split:
        input:
            snp = snp_file,
            bam = mapping_prg+"/{sample}.bam"
        output:
            expand(mapping_prg+"/{{sample}}.{suffix}.bam", suffix = ['allele_flagged', 'genome1', 'genome2', 'unassigned'])
        params:
            paired = SNPparam
        log:
            "allelic_bams/{sample}_SNPsplit.log"
        shell:
            SNPsplit_path + "SNPsplit"
            " {params.paired} --samtools_path "+samtools_path+"samtools"
            " --snp_file {input.snp} {input.bam} 2> {log}"

# move the allele-specific bams to another folder
if mapping_prg != "Bowtie2":
    rule movebams:
        input:
            mapping_prg+"/{sample}.{suffix}.bam"
        output:
            "allelic_bams/{sample}.{suffix}.unsorted.bam"
        shell:
            "mv {input} {output}"
else:
    rule movebams:
        input:
            "filtered_bam/{sample}.filtered.{suffix}.bam"
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
    log: "allelic_bams/{sample}.{suffix}.BAMsorting.log"
    shell:
        samtools_path+"samtools sort -@ {threads} -T ${{TMPDIR}} -O bam -o {output} {input} 2> {log}"

# index the sorted files
rule BAMindex_allelic:
    input:
        "allelic_bams/{sample}.{suffix}.sorted.bam"
    output:
        "allelic_bams/{sample}.{suffix}.sorted.bam.bai"
    shell:
        samtools_path+"samtools index {input}"
