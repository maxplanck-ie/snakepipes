
# SNPsplit rule
if paired:
    SNPparam = '--paired'
else:
    SNPparam = ''

## get input bam depending on the mapping prog (use filtered bam in case of chip-seq data)
if mapping_prg != "Bowtie2":
    SNPbam_files = expand(mapping_prg+"/{sample}.bam", sample = samples)
else:
    SNPbam_files = expand("filtered_bam/{sample}.filtered.bam", sample = samples)

rule snp_split:
    input:
        snp = snp_file,
        bam = SNPbam_files
    output:
        expand(mapping_prg+"/{sample}_{suffix}.bam", sample = samples,
                     suffix = ['allele_flagged', 'genome1', 'genome2', 'unassigned'])
    params:
        paired = SNPparam
    log:
        expand("allelic_bams/{sample}_SNPsplit.log", sample = samples)
    shell:
        SNPsplit_path + "SNPsplit"
        " {params.paired} --samtools_path "+samtools_path+"samtools"
        " --snp_file {input.snp} {input.bam} 2> {log}"

# move the allele-specific bams to another folder
rule movebams:
    input:
        expand(mapping_prg+"/{sample}_{suffix}.bam", sample = samples,
                     suffix = ['allele_flagged', 'genome1', 'genome2', 'unassigned'])
    output:
        temp("allelic_bams/{sample}_{suffix}.bam")
    shell:
        "mv {input} {output}"

rule BAMsort_allelic:
    input: expand("allelic_bams/{sample}_{suffix}.bam", sample = samples, suffix = ['allele_flagged', 'genome1', 'genome2', 'unassigned'])
    output:
        "allelic_bams/{sample}_{suffix}.sorted.bam"
    threads:
        12
    log: "allelic_bams/{sample}.BAMsorting.log"
    shell:
        samtools_path+"samtools sort -@ {threads} -T ${{TMPDIR}} -O bam -o {output} {input} 2> {log}"

rule BAMindex_allelic:
    input:
        expand("allelic_bams/{sample}_{suffix}.sorted.bam", sample = samples,
                     suffix = ['allele_flagged', 'genome1', 'genome2', 'unassigned'])
    output:
        "allelic_bams/{sample}_{suffix}.sorted.bam.bai"
    shell:
        samtools_path+"samtools index {input}"
