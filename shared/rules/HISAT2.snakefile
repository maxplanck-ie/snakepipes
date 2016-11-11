### HISAT2 #####################################################################

if paired:
    rule HISAT2:
        input:
            r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
            r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz"
        output:
            align_summary = "HISAT2/{sample}.HISAT2_summary.txt",
            bam = "HISAT2/{sample}.bam",
            splice = "HISAT2/{sample}.splice_sites.txt",
            met = "HISAT2/{sample}.metrics.txt",
            unconc = "HISAT2/{sample}.un-conc.fastq.gz",
            alconc = "HISAT2/{sample}.al-conc.fastq.gz"
        params:
            hisat_opts = "",
            read_orientation = read_orientation
        benchmark:
            "HISAT2/.benchmark/HISAT2.{sample}.benchmark"
        threads: 8
        shell:
            hisat2_path+"hisat2 "
            "-p {threads} "
            "{params.hisat_opts} {params.read_orientation} "
            "-x "+hisat2_index+" "
            "-1 {input.r1} -2 {input.r2} "
            "--novel-splicesite-outfile {output.splice} "
            "--met-file {output.met} "
            "--un-conc-gz {output.unconc} "
            "--al-conc-gz {output.alconc} "
            "2> {output.align_summary} | "
            ""+samtools_path+"samtools view -Sb - | "
            ""+samtools_path+"samtools sort -m 2G -T ${{TMPDIR}}{wildcards.sample} -@ 2 -O bam - > {output.bam} "
            "&& touch {output.unconc} {output.alconc} "
else:
    pass
    # rule HISAT2:
    #     input:
    #         fastq_dir+"/{sample}.fastq.gz"
    #     output:
    #         align_summary = "Bowtie2/{sample}.Bowtie2_summary.txt",
    #         bam = temp("Bowtie2/{sample}.sorted.bam")
    #     params:
    #         bowtie_opts = ""
    #     benchmark:
    #         "Bowtie2/.benchmark/Bowtie2.{sample}.benchmark"
    #     threads: 24
    #     shell:
    #         bowtie2_path+"bowtie2 "
    #             "-x "+bowtie2_index+" -U {input} "
    #             "{params.bowtie_opts} "
    #             "--rg-id {wildcards.sample} --rg CN:mpi-ie_deep_sequencing_unit "
    #             "--rg DS:{wildcards.sample} --rg PL:ILLUMINA --rg SM:{wildcards.sample} "
    #             "-p {threads} "
    #             "2> {output.align_summary} | "
    #             ""+samtools_path+"samtools view -Sbu - | "
    #             ""+samtools_path+"samtools sort -m 2G -T ${{TMPDIR}}{wildcards.sample} -@ 2 -O bam - > {output.bam}"


### samtools_index #############################################################

rule HISAT2_BAM_index:
    input:
        "HISAT2/{sample}.bam"
    output:
        "HISAT2/{sample}.bam.bai"
    shell:
        samtools_path+"samtools index {input}"
