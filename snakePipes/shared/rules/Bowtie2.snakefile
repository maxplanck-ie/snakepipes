### Bowtie2 ####################################################################
if paired:
    rule Bowtie2:
        input:
            r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
            r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz"
        output:
            align_summary = "Bowtie2/{sample}.Bowtie2_summary.txt",
            bam = temp("Bowtie2/{sample}.sorted.bam")# removing since we keep the sambamba output (dupmarked)
        params:
            bowtie2_index=bowtie2_index,
            bowtie_opts = str(bowtie_opts or ''),
            mateOrientation = mateOrientation,
            insertSizeMax = insertSizeMax
        benchmark:
            "Bowtie2/.benchmark/Bowtie2.{sample}.benchmark"
        threads: 24  # 1G per core
        conda: CONDA_DNA_MAPPING_ENV
        shell: """
            MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
            bowtie2 \
            -X {params.insertSizeMax} \
            -x {params.bowtie2_index} -1 {input.r1} -2 {input.r2} \
            {params.bowtie_opts} {params.mateOrientation} \
            --rg-id {wildcards.sample} --rg CN:mpi-ie_deep_sequencing_unit \
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
            bowtie_opts = str(bowtie_opts or '')
        benchmark:
            "Bowtie2/.benchmark/Bowtie2.{sample}.benchmark"
        threads: 24  # 1G per core
        conda: CONDA_DNA_MAPPING_ENV
        shell: """
            MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
            bowtie2 \
            -x {params.bowtie2_index} -U {input} \
            {params.bowtie_opts} \
            --rg-id {wildcards.sample} --rg CN:mpi-ie_deep_sequencing_unit \
            --rg DS:{wildcards.sample} --rg PL:ILLUMINA --rg SM:{wildcards.sample} \
            -p {threads} \
            2> {output.align_summary} | \
            samtools view -Sbu - | \
            samtools sort -m 2G -T $MYTEMP/{wildcards.sample} -@ 2 -O bam - > {output.bam};
            rm -rf $MYTEMP
            """
