### Bowtie2 ####################################################################

if paired:
    rule Bowtie2:
        input:
            r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
            r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz"
        output:
            align_summary = "Bowtie2/{sample}.Bowtie2_summary.txt",
            bam = temp("Bowtie2/{sample}.sorted.bam")
        params:
            bowtie_opts = "-X 1000",
            mate_orientation = mate_orientation
        benchmark:
            "Bowtie2/.benchmark/Bowtie2.{sample}.benchmark"
        threads: 24
        shell:
            bowtie2_path+"bowtie2 "
            "-x "+bowtie2_index+" -1 {input.r1} -2 {input.r2} "
            "{params.bowtie_opts} {params.mate_orientation} "
            "--rg-id {wildcards.sample} --rg CN:mpi-ie_deep_sequencing_unit --rg DS:{wildcards.sample} --rg PL:ILLUMINA --rg SM:{wildcards.sample} "
            "-p {threads} "
            "2> {output.align_summary} | "
            ""+samtools_path+"samtools view -Sb - | "
            ""+samtools_path+"samtools sort -m 2G -T {wildcards.sample} -@ 2 -O bam - > {output.bam}"
else:
    rule Bowtie2:
        input:
            fastq_dir+"/{sample}.fastq.gz"
        output:
            align_summary = "Bowtie2/{sample}.Bowtie2_summary.txt",
            bam = temp("Bowtie2/{sample}.sorted.bam")
        params:
            bowtie_opts = ""
        benchmark:
            "Bowtie2/.benchmark/Bowtie2.{sample}.benchmark"
        threads: 24
        shell:
            bowtie2_path+"bowtie2 "
                "-x "+bowtie2_index+" -U {input} "
                "{params.bowtie_opts} "
                "--rg-id {wildcards.sample} --rg CN:mpi-ie_deep_sequencing_unit --rg DS:{wildcards.sample} --rg PL:ILLUMINA --rg SM:{wildcards.sample} "
                "-p {threads} "
                "2> {output.align_summary} | "
                ""+samtools_path+"samtools view -Sb - | "
                ""+samtools_path+"samtools sort -m 2G -T {wildcards.sample} -@ 2 -O bam - > {output.bam}"


### Picard MarkDuplicates ######################################################

rule MarkDuplicates:
    input:
        "Bowtie2/{sample}.sorted.bam"
    output:
        bam = "Bowtie2/{sample}.bam",
        txt = "Picard_qc/MarkDuplicates/{sample}.mark_duplicates_metrics.txt"
    log:
        "Picard_qc/logs/MarkDuplicates.{sample}.log"
    benchmark:
        "Picard_qc/.benchmark/MarkDuplicates.{sample}.benchmark"
    threads: 2 # Java performs parallel garbage collection
    shell:
        "java -Xmx8g -jar "+picard_path+"picard.jar MarkDuplicates "
            "MAX_FILE_HANDLES=1000 "
            "INPUT={input} "
            "OUTPUT={output.bam} "
            "METRICS_FILE={output.txt} "
            "ASSUME_SORTED=true "
            "VERBOSITY=WARNING "
            "VALIDATION_STRINGENCY=LENIENT "
            "&> {log} "


### samtools_index #############################################################

rule samtools_index:
    input:
        "Bowtie2/{sample}.bam"
    output:
        "Bowtie2/{sample}.bam.bai"
    shell:
        samtools_path+"samtools index {input}"
