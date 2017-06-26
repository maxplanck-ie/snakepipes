
## SNPFILE still not defined
# SNPsplit rule
if paired:
    SNPparam = '--paired'
else:
    SNPparam = ''

rule snp_split:
    input:
        snp = snp_file,
        sam = mapping_prg+"/{sample}_nsorted.bam"
    output:
        bam = mapping_prg+"/{sample}.allele_flagged.bam",
        allele1 = mapping_prg+"/{sample}.genome1.bam",
        allele2 = mapping_prg+"/{sample}.genome2.bam",
        unassigned = mapping_prg+"/{sample}.unassigned.bam"
    params:
        paired = SNPparam
    log:
        mapping_prg+"/{sample}_SNPsplit.log"
    shell:
        SNPsplit_path + "SNPsplit"
        " {params.paired} --no_sort --samtools_path "+samtools_path+"samtools"
        " --snp_file {input.snp} {input.sam} 2> {log}"

rule bamsort:
    input:
        mapping_prg+"/{sample}_nsorted.bam"
    output:
        mapping_prg+"/{sample}.bam"
    threads:
        12
    log: mapping_prg+"/{sample}.BAMsorting.log"
    shell:
        samtools_path+"samtools sort -@ {threads} -T ${{TMPDIR}} -O bam -o {output} {input} 2> {log}"

rule BAM_index:
    input:
        mapping_prg+"/{sample}.bam"
    output:
        mapping_prg+"/{sample}.bam.bai"
    shell:
        samtools_path+"samtools index {input}"
