rule library_type_downsample:
    input:  indir+"/{sample}{read}.fastq.gz"
    output: "library_type/{sample}{read}.fastq.gz"
    params:
        num = 200000
    benchmark:  "FASTQ/.benchmark/library_type_downsample.{sample}{read}.benchmark"
    threads: 4
    shell:  "bash "+os.path.join(workflow_tools,"fastq_head.sh")+" {threads} {params.num} {input} {output}"


if paired:
    rule library_type_Bowtie2:
        input:
            r1 = "library_type/{sample}"+reads[0]+".fastq.gz",
            r2 = "library_type/{sample}"+reads[1]+".fastq.gz"
        output:
            align_summary = "library_type/{sample}.align_summary.txt",
            accepted_hits = "library_type/{sample}.bam"
        params:
            tmpprefix = "{sample}"
        threads: 4
        benchmark:  "library_type/.benchmark/library_type_Bowtie2.{sample}.benchmark"
        shell:  bowtie2_path+" -p {threads} "
                "-x "+bowtie2_index+" -1 {input.r1} -2 {input.r2} 2> {output.align_summary} | "
                ""+samtools_path+" view -Sb - | "
                ""+samtools_path+" sort -O bam -m 2G -@ 2 -T {params.tmpprefix} - > {output.accepted_hits} "
else:
    rule library_type_Bowtie2:
        input: "library_type/{sample}.fastq.gz"
        output:
            align_summary = "library_type/{sample}.align_summary.txt",
            accepted_hits = "library_type/{sample}.bam"
        params:
            tmpprefix = "{sample}"
        threads: 4
        benchmark:  "library_type/.benchmark/library_type_Bowtie2.{sample}.benchmark"
        shell:  bowtie2_path+" -p {threads} "
                "-x "+bowtie2_index+" -U {input} 2> {output.align_summary} | "
                ""+samtools_path+" view -Sb - | "
                ""+samtools_path+" sort -O bam -m 2G -@ 2 -T {params.tmpprefix} - > {output.accepted_hits} "


rule library_type_RSeQC_infer_experiment:
    input:
        bam = "library_type/{sample}.bam",
        bed = genes_bed
    output: "library_type/{sample}.infer_experiment.txt"
    benchmark:  "library_type/.benchmark/library_type_RSeQC_infer_experiment.{sample}.benchmark"
    log:    "library_type/logs/{sample}.library_type_RSeQC_infer_experiment.log"
    shell:  "bash "+os.path.join(rseqc_path, "activate")+" && "
            ""+os.path.join(rseqc_path,"infer_experiment.py")+" -i {input.bam} -r {input.bed} >{output} 2>{log} "


rule library_type:
    input:  [os.path.join("library_type",x+".infer_experiment.txt") for x in samples]
    output: "library_type/library_type.txt"
    log:    "library_type/logs/library_type.log"
    run:
        input = " ".join(input)
        ##print(input)
        shell( "python "+os.path.join(shared_dir_path,"tools","library_type.py")+" {input} >{output} 2>{log}" )
