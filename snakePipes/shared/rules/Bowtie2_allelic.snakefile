## get basename of the bt2 index
import os

def getbw_idxbase(file):
    base = os.path.basename(file)
    idxbase = re.sub('.[0-9*].bt2','',base)
    fpath = os.path.dirname(file) + "/" + idxbase
    return(fpath)

### Bowtie2 ####################################################################
if aligner == "Bowtie2":
    if pairedEnd:
        rule Bowtie2_allele:
            input:
                r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
                r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz",
                index = bowtie2_index_allelic
            output:
                align_summary = aligner+"/{sample}.Bowtie2_summary.txt",
                bam = temp(aligner+"/{sample}.sorted.bam")
            params:
                alignerOpts = str(alignerOpts or ''),
                mateOrientation = mateOrientation,
                insertSizeMax = insertSizeMax,
                idxbase = getbw_idxbase(bowtie2_index_allelic),
                tempDir = tempDir
            benchmark:
                aligner+"/.benchmark/Bowtie2.{sample}.benchmark"
            threads: lambda wildcards: 24 if 24<max_thread else max_thread  # 1G per core
            conda: CONDA_DNA_MAPPING_ENV
            shell: """
                TMPDIR={params.tempDir}
                MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
                bowtie2 \
                -X {params.insertSizeMax} \
                -x {params.idxbase} -1 {input.r1} -2 {input.r2} \
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
        rule Bowtie2_allele:
            input:
                r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
                index = bowtie2_index_allelic
            output:
                align_summary = aligner+"/{sample}.Bowtie2_summary.txt",
                bam = temp(aligner+"/{sample}.sorted.bam")
            params:
                alignerOpts = str(alignerOpts or ''),
                idxbase = getbw_idxbase(bowtie2_index_allelic),
                tempDir = tempDir
            benchmark:
                aligner+"/.benchmark/Bowtie2.{sample}.benchmark"
            threads: lambda wildcards: 24 if 24<max_thread else max_thread  # 1G per core
            conda: CONDA_DNA_MAPPING_ENV
            shell: """
                TMPDIR={params.tempDir}
                MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
                bowtie2 \
                -x {params.idxbase} -U {input.r1} \
                --reorder \
                {params.alignerOpts} \
                --rg-id {wildcards.sample} \
                --rg DS:{wildcards.sample} --rg PL:ILLUMINA --rg SM:{wildcards.sample} \
                -p {threads} \
                2> {output.align_summary} | \
                samtools view -Sbu - | \
                samtools sort -m 2G -T $MYTEMP/{wildcards.sample} -@ 2 -O bam - > {output.bam};
                rm -rf $MYTEMP
                """
else:
    print("Only bowtie2 implemented")
