if pairedEnd:
    rule bwa:
        input:
            r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
            r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz"
        output:
            align_summary = "bwa/{sample}.bwa_summary.txt", #samtools flagstat
            bam = temp("bwa/{sample}.sorted.bam")
        log: "bwa/logs/{sample}.sort.log"
        params:
            bwa_index = bwa_index,
            alignerOpts = str(alignerOpts or ''),
            tempDir = tempDir
        threads: 24
        conda: CONDA_DNA_MAPPING_ENV 
        shell:"""
            TMPDIR={params.tempDir}
            MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
            bwa mem \
            -t {threads} \
            -R '@RG\\tID:{wildcards.sample}\\tDS:{wildcards.sample}\\tPL:ILLUMINA\\tSM:{wildcards.sample}' \
            {params.bwa_index} {input.r1} {input.r2} | \
            samtools view -Sb - | \
            samtools sort -m 2G -@ 2 -O bam - > {output.bam} 2> {log};
            rm -rf $MYTEMP
            samtools flagstat {output.bam} > {output.align_summary}
        """
else:
    rule bwa:
        input:
            fastq_dir+"/{sample}"+reads[0]+".fastq.gz"
        output:
            align_summary = "bwa/{sample}.bwa_summary.txt", #samtools flagstat
            bam = temp("bwa/{sample}.sorted.bam")
        log: "bwa/logs/{sample}.sort.log"
        params:
            bwa_index = bwa_index,
            alignerOpts = str(alignerOpts or ''),
            tempDir = tempDir
        threads : 24 
        conda: CONDA_DNA_MAPPING_ENV
        shell: """
            TMPDIR={params.tempDir}
            MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
            bwa mem \
            -t {threads} \
            -R '@RG\\tID:{wildcards.sample}\\tDS:{wildcards.sample}\\tPL:ILLUMINA\\tSM:{wildcards.sample}' \
            {params.bwa_index} {input} | \
            samtools view -Sbu - | \
            samtools sort -m 2G -T $MYTEMP/{wildcards.sample} -@ 2 -O bam - > {output.bam} 2> {log};
            rm -rf $MYTEMP
            samtools flagstat {output.bam} > {output.align_summary}
        """
