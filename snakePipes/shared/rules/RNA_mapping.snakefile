def getHISAT_libtype(pairedEnd, libraryType):
    """
    Convert from a featureCounts library type to a HISAT2 option string
    """
    if pairedEnd:
        if libraryType == 1:
            return "--rna-strandness FR"
        elif libraryType == 2:
            return "--rna-strandness RF"
        else:
            return ""
    else:
        if libraryType == 1:
            return "--rna-strandness F"
        elif libraryType == 2:
            return "--rna-strandness R"
        else:
            return ""



### HISAT2 #####################################################################

if aligner.upper().find("HISAT2") >=0:
    if pairedEnd:
        rule HISAT2:
            input:
                r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
                r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz",
            output:
                align_summary = aligner+"/{sample}.HISAT2_summary.txt",
                bam = temp(aligner+"/{sample}.sorted.bam"),
                splice = aligner+"/{sample}/splice_sites.txt",
                met = aligner+"/{sample}/metrics.txt"
            params:
                lib_type = getHISAT_libtype(pairedEnd, libraryType),
                input_splice = known_splicesites,
                alignerOptions = str(alignerOptions or ''),
                samsort_memory = '2G',
                idx = hisat2_index,
                tempDir = tempDir
            benchmark:
                aligner+"/.benchmark/HISAT2.{sample}.benchmark"
            threads: lambda wildcards: 10 if 10<max_thread else max_thread
            conda: CONDA_RNASEQ_ENV
            shell: """
                TMPDIR={params.tempDir}
                MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
                hisat2 -p {threads} {params.alignerOptions} \
                    {params.lib_type} -x {params.idx} \
                    --known-splicesite-infile {params.input_splice} \
                    -1 {input.r1} -2 {input.r2} \
                    --novel-splicesite-outfile {output.splice} \
                    --met-file {output.met} 2> {output.align_summary} \
                | samtools sort -m {params.samsort_memory} -T $MYTEMP/{wildcards.sample} -@ {threads} -O bam -o {output.bam} -;
                rm -rf $MYTEMP
                """
    else:
        rule HISAT2:
            input:
                fastq_dir+"/{sample}"+reads[0]+".fastq.gz"
            output:
                align_summary = aligner+"/{sample}.HISAT2_summary.txt",
                bam = temp(aligner+"/{sample}.sorted.bam"),
                splice = aligner+"/{sample}/splice_sites.txt",
                met = aligner+"/{sample}/metrics.txt"
            params:
                lib_type = getHISAT_libtype(pairedEnd, libraryType),
                input_splice = known_splicesites,
                alignerOptions = str(alignerOptions or ''),
                samsort_memory = '2G',
                idx = hisat2_index,
                tempDir = tempDir
            benchmark:
                aligner+"/.benchmark/HISAT2.{sample}.benchmark"
            threads: lambda wildcards: 10 if 10<max_thread else max_thread
            conda: CONDA_RNASEQ_ENV
            shell: """
                TMPDIR={params.tempDir}
                MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX)
                hisat2 -p {threads} {params.alignerOptions} \
                    {params.lib_type} -x {params.idx} \
                    --known-splicesite-infile {params.input_splice} \
                    -U {input[0]} \
                    --novel-splicesite-outfile {output.splice} \
                    --met-file {output.met} 2> {output.align_summary} \
                | samtools sort -m {params.samsort_memory} -T $MYTEMP/{wildcards.sample} -@ {threads} -O bam -o {output.bam} -
                rm -rf $MYTEMP
                """
elif aligner.upper().find("STAR") >=0:
    if pairedEnd:
        rule STAR:
            input:
                r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
                r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz"
            output:
                bam = temp(aligner+"/{sample}.sorted.bam")
            params:
                alignerOptions = str(alignerOptions or ''),
                gtf = genes_gtf,
                index = star_index,
                prefix = aligner+"/{sample}/{sample}.",
                samsort_memory = '2G',
                sample_dir = aligner+"/{sample}",
                samtools_threads = 5,
                tempDir = tempDir
            benchmark:
                aligner+"/.benchmark/STAR.{sample}.benchmark"
            threads: lambda wildcards: 20 if 20<max_thread else max_thread  # 3.2G per core
            conda: CONDA_RNASEQ_ENV
            shell: """
                TMPDIR={params.tempDir}
                MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX)
                ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
                STAR --runThreadN {threads} \
                    {params.alignerOptions} \
                    --sjdbOverhang 100 \
                    --outSAMunmapped Within \
                    --outSAMtype BAM Unsorted \
                    --outStd BAM_Unsorted \
                    --sjdbGTFfile {params.gtf} \
                    --genomeDir {params.index} \
                    --readFilesIn  {input.r1} {input.r2} \
                    --readFilesCommand 'gunzip -c' \
                    --outFileNamePrefix {params.prefix} \
                | samtools sort -m {params.samsort_memory} -T $MYTEMP/{wildcards.sample} -@ {params.samtools_threads} -O bam -o {output.bam} -
                rm -rf $MYTEMP
                """
    else:
        rule STAR:
            input:
                fastq_dir+"/{sample}"+reads[0]+".fastq.gz"
            output:
                bam = temp(aligner+"/{sample}.sorted.bam")
            params:
                alignerOptions = str(alignerOptions or ''),
                gtf = genes_gtf,
                index = star_index,
                prefix = aligner+"/{sample}/{sample}.",
                samsort_memory = '2G',
                sample_dir = aligner+"/{sample}",
                samtools_threads = 5,
                tempDir = tempDir
            benchmark:
                aligner+"/.benchmark/STAR.{sample}.benchmark"
            threads: lambda wildcards: 20 if 20<max_thread else max_thread  # 3.2G per core
            conda: CONDA_RNASEQ_ENV
            shell: """
                TMPDIR={params.tempDir}
                MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX)
                ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
                STAR --runThreadN {threads} \
                    {params.alignerOptions} \
                    --sjdbOverhang 100 \
                    --outSAMunmapped Within \
                    --outSAMtype BAM Unsorted \
                    --outStd BAM_Unsorted \
                    --sjdbGTFfile {params.gtf} \
                    --genomeDir {params.index} \
                    --readFilesIn {input} \
                    --readFilesCommand 'gunzip -c' \
                    --outFileNamePrefix {params.prefix} \
                | samtools sort -m {params.samsort_memory} -T $MYTEMP/{wildcards.sample} -@ {params.samtools_threads} -O bam -o {output.bam} -
                rm -rf $MYTEMP
                """
