if paired:
    rule Bowtie2:
        input:
            r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
            r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz"
        output:
            align_summary = "Bowtie2/{sample}.align_summary.txt",
            accepted_hits = "Bowtie2/{sample}.accepted_hits.bam"
        params:
            tmpprefix = "{sample}",
            bowtie_opts = "-X 1000",
            mate_orientation = mate_orientation
        threads: 18
        benchmark:  "Bowtie2/.benchmark/Bowtie2.{sample}.benchmark"
        shell:
            bowtie2_path+" -p {threads} {params.bowtie_opts} {params.mate_orientation} "
            "--rg DS:{params.tmpprefix} --rg SM:{params.tmpprefix} --rg-id mpi-ie --rg CN:deep_sequencing_unit --rg PL:illumina "
            "-x "+Bowtie2Index+" -1 {input.r1} -2 {input.r2} 2> {output.align_summary} | "
            ""+samtools_path+" view -Sb - | "
            ""+samtools_path+" sort -O bam -m 2G -@ 2 -T {params.tmpprefix} - > {output.accepted_hits} "
else:
    rule Bowtie2:
        input: fastq_dir+"/{sample}.fastq.gz"
        output:
            align_summary = "Bowtie2/{sample}.align_summary.txt",
            accepted_hits = "Bowtie2/{sample}.accepted_hits.bam"
        params:
            tmpprefix = "{sample}",
            bowtie_opts = ""
        threads: 18
        benchmark:  "Bowtie2/.benchmark/Bowtie2.{sample}.benchmark"
        shell:  bowtie2_path+" -p {threads} {params.bowtie_opts} "
                "--rg DS:{params.tmpprefix} --rg SM:{params.tmpprefix} --rg-id mpi-ie --rg CN:deep_sequencing_unit --rg PL:illumina "
                "-x "+Bowtie2Index+" -U {input} 2> {output.align_summary} | "
                ""+samtools_path+" view -Sb - | "
                ""+samtools_path+" sort -O bam -m 2G -@ 2 -T {params.tmpprefix} - > {output.accepted_hits} "


rule MarkDuplicates:
    input:  "Bowtie2/{sample}.accepted_hits.bam"
    output:
        bam = "Bowtie2/{sample}.bam",
        dat = "Bowtie2/{sample}.MarkDuplicates.dat"
    params:
        sample = "{sample}"
    log:    "Bowtie2/log/{sample}.MarkDuplicates.log"
    benchmark:  "Bowtie2/.benchmark/Picard_MarkDuplicates.{sample}.benchmark"
    shell:  "java -Xmx4g -jar "+picard_path+" MarkDuplicates "
            "I={input} O={output.bam} M={output.dat} "
            "VALIDATION_STRINGENCY=LENIENT VERBOSITY=WARNING ASSUME_SORTED=true MAX_FILE_HANDLES=1000 "
            "&& rm {input} "
            "2>&1 | tee {log} "


rule MarkDuplicates_bam_index:
    input:  "Bowtie2/{sample}.bam"
    output: "Bowtie2/{sample}.bam.bai"
    shell:  samtools_path+" index {input}"
