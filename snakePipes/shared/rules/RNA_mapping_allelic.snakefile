## allelic mapping using STAR
if config['aligner'] == "STAR":
    if config['pairedEnd']:
        rule STAR_allele:
            input:
                r1 = config['fastq_dir']+"/{sample}"+config['reads'][0]+".fastq.gz",
                r2 = config['fastq_dir']+"/{sample}"+config['reads'][1]+".fastq.gz",
                index = config['star_index_allelic']
            output:
                temp(config['aligner']+"/{sample}.sorted.bam")
            log: config['aligner']+"/logs/{sample}.sort.log"
            params:
                alignerOptions = str(config['alignerOptions'] or ''),
                gtf = config['genes_gtf'],
                prefix = config['aligner']+"/{sample}/{sample}.",
                samsort_memory = '2G',
                idx = os.path.dirname(config['star_index_allelic']),
                sample_dir = config['aligner']+"/{sample}",
                tempDir = config['tempDir']
            benchmark:
                config['aligner']+"/.benchmark/STAR.{sample}.benchmark"
            threads: lambda wildcards: 12 if 12<config['max_thread'] else config['max_thread']
            conda: config['CONDA_RNASEQ_ENV']
            shell: """
                TMPDIR={params.tempDir}
                MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
                ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
                STAR {params.alignerOptions} \
                    --runThreadN {threads} \
                    --genomeDir {params.idx} \
                    --sjdbGTFfile {params.gtf} \
                    --sjdbOverhang 100 \
                    --readFilesIn <(gunzip -c {input.r1}) <(gunzip -c {input.r2}) \
                    --outFileNamePrefix {params.prefix} \
                    --outSAMunmapped Within \
                    --alignEndsType EndToEnd \
                    --outSAMattributes NH HI NM MD \
                    --outSAMtype BAM Unsorted \
                    --outStd BAM_Unsorted \
                    --outFilterMultimapNmax 20 \
                    --alignSJoverhangMin 8 \
                    --alignSJDBoverhangMin 1 \
                    --outFilterMismatchNmax 999 \
                    --alignIntronMin 1 \
                    --alignIntronMax 1000000 \
                    --alignMatesGapMax 1000000 \
                | samtools sort -m {params.samsort_memory} -T $MYTEMP/{wildcards.sample} -@ {threads} -O bam -o {output} - 2> {log};
                rm -rf $MYTEMP
                """
    else:
        rule STAR_allele:
            input:
                r1 = config['fastq_dir']+"/{sample}.fastq.gz",
                index = config['star_index_allelic']
            output:
                temp(config['aligner']+"/{sample}.sorted.bam")
            log: config['aligner']+"/logs/{sample}.sort.log"
            params:
                alignerOptions = str(config['alignerOptions'] or ''),
                gtf = config['genes_gtf'],
                idx = os.path.dirname(config['star_index_allelic']),
                prefix = config['aligner']+"/{sample}/{sample}.",
                samsort_memory = '2G',
                sample_dir = config['aligner']+"/{sample}",
                tempDir = config['tempDir']
            benchmark:
                config['aligner']+"/.benchmark/STAR.{sample}.benchmark"
            threads: lambda wildcards: 12 if 12<config['max_thread'] else config['max_thread']
            conda: config['CONDA_RNASEQ_ENV']
            shell: """
                TMPDIR={params.tempDir}
                MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
                ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
                STAR {params.alignerOptions} \
                    --runThreadN {threads} \
                    --genomeDir {params.idx} \
                    --sjdbGTFfile {params.gtf} \
                    --sjdbOverhang 100 \
                    --readFilesIn <(gunzip -c {input}) \
                    --outFileNamePrefix {params.prefix} \
                    --outSAMunmapped Within \
                    --alignEndsType EndToEnd \
                    --outSAMattributes NH HI NM MD  \
                    --outSAMtype BAM Unsorted \
                    --outStd BAM_Unsorted \
                    --outFilterMultimapNmax 20 \
                    --alignSJoverhangMin 8 \
                    --alignSJDBoverhangMin 1 \
                    --outFilterMismatchNmax 999 \
                    --alignIntronMin 1 \
                    --alignIntronMax 1000000 \
                    --alignMatesGapMax 1000000 \
                | samtools sort -m {params.samsort_memory} -T $MYTEMP/{wildcards.sample} -@ {threads} -O bam -o {output} - 2> {log};
                rm -rf $MYTEMP
                """
else:
    print("Only STAR is implemented for Allele-specific mapping")
