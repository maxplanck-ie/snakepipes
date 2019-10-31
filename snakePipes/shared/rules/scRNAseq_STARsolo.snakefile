##STARsolo

rule STARsolo:
    input:
        r1="originalFASTQ/{sample}"+reads[0]+".fastq.gz",
        r2="originalFASTQ/{sample}"+reads[1]+".fastq.gz"
    output:
        bam = "STARsolo/{sample}.sorted.bam"
    params:
        alignerOptions = str(alignerOptions or ''),
        gtf = genes_gtf,
        index = star_index,
        prefix = "STARsolo/{sample}/{sample}.",
        samsort_memory = '2G',
        sample_dir = "STARsolo/{sample}"
    benchmark:
        aligner+"/.benchmark/STARsolo.{sample}.benchmark"
    threads: 20  # 3.2G per core
    conda: CONDA_scRNASEQ_ENV
    shell: """
        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
        ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
        STAR --runThreadN {threads} \
            {params.alignerOptions} \
            --sjdbOverhang 100 \
            --readFilesCommand zcat \
            --outSAMunmapped Within \
            --outSAMtype BAM Unsorted \
            --outStd BAM_Unsorted \
            --sjdbGTFfile {params.gtf} \
            --genomeDir {params.index} \
            --readFilesIn {input.r2} {input.r1} \
            --outFileNamePrefix {params.prefix} \
	    --soloType CB_UMI_Simple \
	    --soloUMIstart 1 \
	    --soloUMIlen 7 \
	    --soloCBstart 8 \
	    --soloCBlen 7 \
	    --soloCBwhitelist /data/processing/bioinfo-core/celseq_barcodes.384.1col.txt \
	    --soloBarcodeReadLength 0 \
	    --soloCBmatchWLtype Exact \
	    --soloStrand Reverse\
	    --soloUMIdedup Exact \
        | samtools sort -m {params.samsort_memory} -T $MYTEMP/{wildcards.sample} -@ {threads} -O bam -o {output.bam} -;
        rm -rf $MYTEMP
        """

#rule velocyto:
#    input:
#    output:
#    params:
#    threads:
#    conda:
#    shell:
