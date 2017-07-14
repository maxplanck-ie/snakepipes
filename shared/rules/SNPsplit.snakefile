
# SNPsplit rule
if paired:
    SNPparam = '--paired'
else:
    SNPparam = ''

## get input bam depending on the mapping prog (use filtered bam in case of chip-seq data)
if mapping_prg != "Bowtie2":
    SNPbam_files = mapping_prg+"/{sample}.bam"
else:
    SNPbam_files = "filtered_bam/{sample}.filtered.bam"

## define SNP splitting separately for single and dual hybrid genomes??
rule snp_split:
    input:
        snp = snp_file,
        bam = SNPbam_files
    output:
         expand(mapping_prg+"/{{sample}}_{suffix}.bam", suffix = ['allele_flagged', 'genome1', 'genome2', 'unassigned'])
    params:
        paired = SNPparam
    log:
        "allelic_bams/{sample}_SNPsplit.log"
    shell:
        SNPsplit_path + "SNPsplit"
        " {params.paired} --samtools_path "+samtools_path+"samtools"
        " --snp_file {input.snp} {input.bam} 2> {log}"

# move the allele-specific bams to another folder
rule movebams:
    input:
        mapping_prg+"/{sample}_{suffix}.bam"
    output:
        temp("allelic_bams/{sample}_{suffix}.bam")
    shell:
        "mv "+mapping_prg+"{input} {output}"

# sort them
rule BAMsort_allelic:
    input: "allelic_bams/{sample}_{suffix}.bam"
    output:
        "allelic_bams/{sample}_{suffix}.sorted.bam"
    threads:
        12
    log: "allelic_bams/{sample}_{suffix}.BAMsorting.log"
    shell:
        samtools_path+"samtools sort -@ {threads} -T ${{TMPDIR}} -O bam -o {output} {input} 2> {log}"

# index the sorted files
rule BAMindex_allelic:
    input:
        "allelic_bams/{sample}_{suffix}.sorted.bam"
    output:
        "allelic_bams/{sample}_{suffix}.sorted.bam.bai"
    shell:
        samtools_path+"samtools index {input}"
