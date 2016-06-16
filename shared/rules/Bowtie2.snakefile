rule Bowtie2:
    input:
        r1 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
        r2 = "FASTQ/{sample}"+reads[1]+".fastq.gz"
    output:
        align_summary = "Bowtie2/{sample}.align_summary.txt",
        accepted_hits = "Bowtie2/{sample}.bam"
    params:
        tmpprefix = "{sample}",
        bowtie_opts = "-X 1000"
    threads: 18
    benchmark:  "Bowtie2/.benchmark/Bowtie2.{sample}.benchmark"
    shell:  bowtie2_path+" -p {threads} {params.bowtie_opts} "
            "--rg DS:{params.tmpprefix} --rg SM:{params.tmpprefix} --rg-id mpi-ie --rg CN:deep_sequencing_unit --rg PL:illumina "
            "-x "+bowtie2_index+" -1 {input.r1} -2 {input.r2} 2> {output.align_summary} | "
            ""+samtools_path+" view -Sb - | "
            ""+samtools_path+" sort -O bam -m 2G -@ 2 -T {params.tmpprefix} - > {output.accepted_hits} "


rule Bowtie2_bam_index:
    input:  "Bowtie2/{sample}.bam"
    output: "Bowtie2/{sample}.bam.bai"
    shell:  ""+samtools_path+" index {input}"
