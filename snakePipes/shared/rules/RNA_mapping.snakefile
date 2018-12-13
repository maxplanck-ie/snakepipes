def getHISAT_libtype(paired, library_type):
    """
    Convert from a featureCounts library type to a HISAT2 option string
    """
    if paired:
        if library_type == 1:
            return "--rna-strandness FR"
        elif library_type == 2:
            return "--rna-strandness RF"
        else:
            return ""
    else:
        if library_type == 1:
            return "--rna-strandness F"
        elif library_type == 2:
            return "--rna-strandness R"
        else:
            return ""



### HISAT2 #####################################################################

if mapping_prg.upper().find("HISAT2") >=0:
    if paired:
        rule HISAT2:
            input:
                r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
                r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz",
            output:
                align_summary = mapping_prg+"/{sample}.HISAT2_summary.txt",
                bam = temp(mapping_prg+"/{sample}.sorted.bam"),
                splice = mapping_prg+"/{sample}/splice_sites.txt",
                met = mapping_prg+"/{sample}/metrics.txt"
            params:
                lib_type = getHISAT_libtype(paired, library_type),
                input_splice = known_splicesites,
                hisat_options = str(hisat_options or ''),
                samsort_memory = '2G',
                idx = hisat2_index
            benchmark:
                mapping_prg+"/.benchmark/HISAT2.{sample}.benchmark"
            threads: 10
            conda: CONDA_RNASEQ_ENV
            shell: """
                hisat2 -p {threads} {params.hisat_options} \
                    {params.lib_type} -x {params.idx} \
                    --known-splicesite-infile {params.input_splice} \
                    -1 {input.r1} -2 {input.r2} \
                    --novel-splicesite-outfile {output.splice} \
                    --met-file {output.met} 2> {output.align_summary} \
                | samtools sort -m {params.samsort_memory} -T ${{TMPDIR}}{wildcards.sample} -@ {threads} -O bam -o {output.bam} -
                """
    else:
        rule HISAT2:
            input:
                fastq_dir+"/{sample}"+reads[0]+".fastq.gz"
            output:
                align_summary = mapping_prg+"/{sample}.HISAT2_summary.txt",
                bam = temp(mapping_prg+"/{sample}.sorted.bam"),
                splice = mapping_prg+"/{sample}/splice_sites.txt",
                met = mapping_prg+"/{sample}/metrics.txt"
            params:
                lib_type = getHISAT_libtype(paired, library_type),
                input_splice = known_splicesites,
                hisat_options = str(hisat_options or ''),
                samsort_memory = '2G',
                idx = hisat2_index
            benchmark:
                mapping_prg+"/.benchmark/HISAT2.{sample}.benchmark"
            threads: 10
            conda: CONDA_RNASEQ_ENV
            shell: """
                hisat2 -p {threads} {params.hisat_options} \
                    {params.lib_type} -x {params.idx} \
                    --known-splicesite-infile {params.input_splice} \
                    -U {input[0]} \
                    --novel-splicesite-outfile {output.splice} \
                    --met-file {output.met} 2> {output.align_summary} \
                | samtools sort -m {params.samsort_memory} -T ${{TMPDIR}}{wildcards.sample} -@ {threads} -O bam -o {output.bam} -
                """
elif mapping_prg.upper().find("STAR") >=0:
    if paired:
        rule STAR:
            input:
                r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
                r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz"
            output:
                bam = temp(mapping_prg+"/{sample}.sorted.bam")
            params:
                star_options = str(star_options or ''),
                gtf = genes_gtf,
                index = star_index,
                prefix = mapping_prg+"/{sample}/{sample}.",
                samsort_memory = '2G',
                sample_dir = mapping_prg+"/{sample}"
            benchmark:
                mapping_prg+"/.benchmark/STAR.{sample}.benchmark"
            threads: 20  # 3.2G per core
            conda: CONDA_RNASEQ_ENV
            shell: """
                ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
                STAR --runThreadN {threads} \
                    {params.star_options} \
                    --sjdbOverhang 100 \
                    --readFilesCommand zcat \
                    --outSAMunmapped Within \
                    --outSAMtype BAM Unsorted \
                    --outStd BAM_Unsorted \
                    --sjdbGTFfile {params.gtf} \
                    --genomeDir {params.index} \
                    --readFilesIn {input.r1} {input.r2} \
                    --outFileNamePrefix {params.prefix} \
                | samtools sort -m {params.samsort_memory} -T ${{TMPDIR}}{wildcards.sample} -@ {threads} -O bam -o {output.bam} -
                """
    else:
        rule STAR:
            input:
                fastq_dir+"/{sample}"+reads[0]+".fastq.gz"
            output:
                bam = temp(mapping_prg+"/{sample}.sorted.bam")
            params:
                star_options = str(star_options or ''),
                gtf = genes_gtf,
                index = star_index,
                prefix = mapping_prg+"/{sample}/{sample}.",
                samsort_memory = '2G',
                sample_dir = mapping_prg+"/{sample}"
            benchmark:
                mapping_prg+"/.benchmark/STAR.{sample}.benchmark"
            threads: 20  # 3.2G per core
            conda: CONDA_RNASEQ_ENV
            shell: """
                ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
                STAR --runThreadN {threads} \
                    {params.star_options} \
                    --sjdbOverhang 100 \
                    --readFilesCommand zcat \
                    --outSAMunmapped Within \
                    --outSAMtype BAM Unsorted \
                    --outStd BAM_Unsorted \
                    --sjdbGTFfile {params.gtf} \
                    --genomeDir {params.index} \
                    --readFilesIn {input} \
                    --outFileNamePrefix {params.prefix} \
                | samtools sort -m {params.samsort_memory} -T ${{TMPDIR}}{wildcards.sample} -@ {threads} -O bam -o {output.bam} -
                """
