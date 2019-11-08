##STARsolo
##remember that reads are swapped in internals.snakefile!!
rule STARsolo:
    input:
        r1="originalFASTQ/{sample}"+reads[0]+".fastq.gz",
        r2="originalFASTQ/{sample}"+reads[1]+".fastq.gz"
    output:
        bam = "STARsolo/{sample}.sorted.bam"
    params:
        alignerOptions = str(alignerOptions or ''),
        gtf = outdir+"/Annotation/genes.filtered.gtf",
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
            --readFilesIn {input.r1} {input.r2} \
            --outFileNamePrefix {params.prefix} \
	    --soloType CB_UMI_Simple \
	    --soloUMIstart 1 \
	    --soloUMIlen 7 \
	    --soloCBstart 8 \
	    --soloCBlen 7 \
	    --soloCBwhitelist /data/processing/bioinfo-core/celseq_barcodes.384.1col.txt \
	    --soloBarcodeReadLength 0 \
	    --soloCBmatchWLtype Exact \
	    --soloStrand Forward\
	    --soloUMIdedup Exact \
        | samtools sort -m {params.samsort_memory} -T $MYTEMP/{wildcards.sample} -@ {threads} -O bam -o {output.bam} -;
        rm -rf $MYTEMP
        """

rule filter_reads:
    input:
        bamfile = aligner+"/{sample}.sorted.bam",
        bami = aligner+"/{sample}.sorted.bam.bai"
    output:
        bamfile = "filtered_bam/{sample}.filtered.bam",
        bami = "filtered_bam/{sample}.filtered.bam.bai"
    shell: """
           pwd
           ln -s -r {input.bamfile} {output.bamfile} ;
           ln -s -r {input.bami} {output.bami}
           """

#rule velocyto:
#    input:
#    output:
#    params:
#    threads:
#    conda:
#    shell:
