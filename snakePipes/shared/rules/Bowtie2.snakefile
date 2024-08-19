### Bowtie2 ####################################################################
if pairedEnd:
    rule Bowtie2:
        input:
            r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
            r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz"
        output:
            align_summary = "Bowtie2/{sample}.Bowtie2_summary.txt",
            bam = temp("Bowtie2/{sample}.sorted.bam")# removing since we keep the sambamba output (dupmarked)
        params:
            bowtie2_index=bowtie2_index,
            alignerOpts = str(alignerOpts or ' ') if not cutntag else "  --local --very-sensitive-local "\
            "--no-mixed --no-discordant --phred33 -I 10 ",
            mateOrientation = mateOrientation,
            insertSizeMax = str(insertSizeMax or ' ') if not cutntag else " 700 ",
            tempDir = tempDir
        benchmark:
            "Bowtie2/.benchmark/Bowtie2.{sample}.benchmark"
        threads: lambda wildcards: 24 if 24<max_thread else max_thread  # 1G per core
        conda: CONDA_DNA_MAPPING_ENV
        shell: """
            TMPDIR={params.tempDir}
            MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
            bowtie2 \
            -X {params.insertSizeMax} \
            -x {params.bowtie2_index} -1 {input.r1} -2 {input.r2} \
            {params.alignerOpts} {params.mateOrientation} \
            --rg-id {wildcards.sample} \
            --rg DS:{wildcards.sample} --rg PL:ILLUMINA --rg SM:{wildcards.sample} \
            -p {threads} \
            2> {output.align_summary} | \
            samtools view -Sb - | \
            samtools sort -m 2G -T $MYTEMP/{wildcards.sample} -@ 2 -O bam - > {output.bam};
            rm -rf $MYTEMP
            """
else:
    rule Bowtie2:
        input:
            fastq_dir+"/{sample}"+reads[0]+".fastq.gz"
        output:
            align_summary = "Bowtie2/{sample}.Bowtie2_summary.txt",
            bam = temp("Bowtie2/{sample}.sorted.bam")
        params:
            bowtie2_index=bowtie2_index,
            alignerOpts = str(alignerOpts or ''),
            tempDir = tempDir
        benchmark:
            "Bowtie2/.benchmark/Bowtie2.{sample}.benchmark"
        threads: lambda wildcards: 24 if 24<max_thread else max_thread  # 1G per core
        conda: CONDA_DNA_MAPPING_ENV
        shell: """
            TMPDIR={params.tempDir}
            MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
            bowtie2 \
            -x {params.bowtie2_index} -U {input} \
            {params.alignerOpts} \
            --rg-id {wildcards.sample} \
            --rg DS:{wildcards.sample} --rg PL:ILLUMINA --rg SM:{wildcards.sample} \
            -p {threads} \
            2> {output.align_summary} | \
            samtools view -Sbu - | \
            samtools sort -m 2G -T $MYTEMP/{wildcards.sample} -@ 2 -O bam - > {output.bam};
            rm -rf $MYTEMP
            """
