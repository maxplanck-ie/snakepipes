


## allelic mapping using STAR
if mapping_prg == "STAR":
    if paired:
        rule STAR_allele:
            input:
                r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
                r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz",
                index = star_index_allelic
            output:
                mapping_prg+"/{sample}.bam"
            params:
                star_options = str(star_options or ''),
                gtf = genes_gtf,
                prefix = mapping_prg+"/{sample}/{sample}.",
                idx = os.path.dirname(star_index_allelic),
                sample_dir = mapping_prg+"/{sample}"
            benchmark:
                mapping_prg+"/.benchmark/STAR.{sample}.benchmark"
            threads: 12
            conda: CONDA_RNASEQ_ENV
            shell:
                " ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} ) && "
                " STAR"
                " {params.star_options}"
                " --runThreadN {threads}"
                " --genomeDir {params.idx}"
                " --sjdbGTFfile {params.gtf}"
                " --sjdbOverhang 100"
                " --readFilesCommand zcat"
                " --readFilesIn {input.r1} {input.r2}"
                # By default, this prefix is "./".
                " --outFileNamePrefix {params.prefix}"
                " --outSAMunmapped Within"
                # Recommended settings for SNPsplit
                " --alignEndsType EndToEnd"
                " --outSAMattributes NH HI NM MD"
                " --outSAMtype BAM SortedByCoordinate"
                # Additional params added
                " --outFilterMultimapNmax 20"
                " --alignSJoverhangMin 8"
                " --alignSJDBoverhangMin 1"
                " --outFilterMismatchNmax 999"
                " --alignIntronMin 1"
                " --alignIntronMax 1000000"
                " --alignMatesGapMax 1000000"
                " && mv {params.prefix}Aligned.sortedByCoord.out.bam {output} "
    else:
        rule STAR_allele:
            input:
                r1 = fastq_dir+"/{sample}.fastq.gz",
                index = star_index_allelic
            output:
                mapping_prg+"/{sample}.bam"
            params:
                star_options = str(star_options or ''),
                gtf = genes_gtf,
                idx = os.path.dirname(star_index_allelic),
                prefix = mapping_prg+"/{sample}/{sample}.",
                sample_dir = mapping_prg+"/{sample}"
            benchmark:
                mapping_prg+"/.benchmark/STAR.{sample}.benchmark"
            threads: 12
            conda: CONDA_RNASEQ_ENV
            shell:
                " ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} ) && "
                " STAR"
                " {params.star_options}"
                " --runThreadN {threads}"
                " --genomeDir {params.idx}"
                " --sjdbGTFfile {params.gtf}"
                " --sjdbOverhang 100"
                " --readFilesCommand zcat"
                " --readFilesIn {input}"
                # By default, this prefix is "./".
                " --outFileNamePrefix {params.prefix}"
                " --outSAMunmapped Within"
                # Recommended settings for SNPsplit
                " --alignEndsType EndToEnd"
                " --outSAMattributes NH HI NM MD"
                " --outSAMtype BAM SortedByCoordinate"
                # Additional params added
                " --outFilterMultimapNmax 20"
                " --alignSJoverhangMin 8"
                " --alignSJDBoverhangMin 1"
                " --outFilterMismatchNmax 999"
                " --alignIntronMin 1"
                " --alignIntronMax 1000000"
                " --alignMatesGapMax 1000000"
                " && mv {params.prefix}Aligned.sortedByCoord.out.bam {output}"
else:
    print("Only STAR is implemented for Allele-specific mapping")


#### INDEX the mapped files
rule BAM_index:
    input:
        mapping_prg+"/{sample}.bam"
    output:
        mapping_prg+"/{sample}.bam.bai"
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"
