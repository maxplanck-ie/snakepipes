## get basename of the bt2 index
import os

def getbw_idxbase(file):
    base = os.path.basename(file)
    idxbase = re.sub('.[0-9*].bt2','',base)
    fpath = os.path.dirname(file) + "/" + idxbase
    return(fpath)

### Bowtie2 ####################################################################
if mapping_prg == "Bowtie2":
    if paired:
        rule Bowtie2_allele:
            input:
                r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
                r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz",
                index = bowtie2_index_allelic
            output:
                align_summary = mapping_prg+"/{sample}.Bowtie2_summary.txt",
                bam = temp(mapping_prg+"/{sample}.sorted.bam")
            params:
                bowtie_opts = str(bowtie_opts or ''),
                mate_orientation = mate_orientation,
                insert_size_max = insert_size_max,
                idxbase = getbw_idxbase(bowtie2_index_allelic)
            benchmark:
                mapping_prg+"/.benchmark/Bowtie2.{sample}.benchmark"
            threads: 24
            conda: CONDA_DNA_MAPPING_ENV
            shell:
                "bowtie2"
                " -X {params.insert_size_max}"
                " -x {params.idxbase} -1 {input.r1} -2 {input.r2}"
                " {params.bowtie_opts} {params.mate_orientation}"
                " --rg-id {wildcards.sample} --rg CN:mpi-ie_deep_sequencing_unit"
                " --rg DS:{wildcards.sample} --rg PL:ILLUMINA --rg SM:{wildcards.sample}"
                " -p {threads}"
                " 2> {output.align_summary} |"
                "samtools view -Sb - |"
                "samtools sort -m 2G -T ${{TMPDIR}}{wildcards.sample} -@ 2 -O bam - > {output.bam}"
    else:
        rule Bowtie2_allele:
            input:
                r1 = fastq_dir+"/{sample}.fastq.gz",
                index = bowtie2_index_allelic
            output:
                align_summary = mapping_prg+"/{sample}.Bowtie2_summary.txt",
                bam = temp(mapping_prg+"/{sample}.sorted.bam")
            params:
                bowtie_opts = str(bowtie_opts or ''),
                idxbase = getbw_idxbase(bowtie2_index_allelic)
            benchmark:
                mapping_prg+"/.benchmark/Bowtie2.{sample}.benchmark"
            threads: 24
            conda: CONDA_DNA_MAPPING_ENV
            shell:
                "bowtie2"
                " -x {params.idxbase} -U {input.r1}"
                " --reorder"
                " {params.bowtie_opts}"
                " --rg-id {wildcards.sample} --rg CN:mpi-ie_deep_sequencing_unit "
                " --rg DS:{wildcards.sample} --rg PL:ILLUMINA --rg SM:{wildcards.sample} "
                " -p {threads}"
                " 2> {output.align_summary} |"
                "samtools view -Sbu - |"
                "samtools sort -m 2G -T ${{TMPDIR}}{wildcards.sample} -@ 2 -O bam - > {output.bam}"
else:
    print("Only bowtie2 implemented")
