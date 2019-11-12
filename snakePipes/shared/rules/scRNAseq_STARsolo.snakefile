##STARsolo
##remember that reads are swapped in internals.snakefile!!
rule STARsolo:
    input:
        r1="originalFASTQ/{sample}"+reads[0]+".fastq.gz",
        r2="originalFASTQ/{sample}"+reads[1]+".fastq.gz",
        annot="Annotation/genes.filtered.gtf"
    output:
        bam = "STARsolo/{sample}/{sample}.star.bam"
    params:
        alignerOptions = str(alignerOptions or ''),
        gtf = outdir+"/Annotation/genes.filtered.gtf",
        index = star_index,
        prefix = "STARsolo/{sample}/{sample}.star",
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
	    --soloUMIdedup Exact
 
            rm -rf $MYTEMP
         """


#rule sort_bam:
#    input:
#        bamfile="STARsolo/{sample}/{sample}.star.bam"
#    output:
#        bamfile = aligner+"/{sample}.sorted.bam",
#        bami = aligner+"/{sample}.sorted.bam.bai"
#    conda: CONDA_scRNASEQ_ENV
#    shell: """
#        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX)
#        samtools sort -m 2g -T $MYTEMP/{wildcards.sample} -@ 4 -O bam -o {output.bamfile} {input.bamfile}
#        samtools index {output.bamfile}
#        rm -rf $MYTEMP
#    """


#rule filter_bam:
#    input:
#        bamfile = aligner+"/{sample}.sorted.bam",
#        bami = aligner+"/{sample}.sorted.bam.bai"
#    output:
#        bamfile = "filtered_bam/{sample}.filtered.bam",
#        bami = "filtered_bam/{sample}.filtered.bam.bai"
#    shell: """
#           pwd
#           ln -s -r {input.bamfile} {output.bamfile} ;
#           ln -s -r {input.bami} {output.bami}
#           """

#the barcode whitelist is currently passed in although it's not tested if it's actually necessery as it was already provided to STARsolo
#gtf mask is not used as a filtered gtf is passed in
#no metadata table is provided

#rule velocyto:
#    input:
#        bc = "/data/processing/bioinfo-core/celseq_barcodes.384.1col.txt",
#        gtf = outdir+"/Annotation/genes.filtered.gtf",
#        bam = expand("STARsolo/{sample}.sorted.bam",sample=samples)
#    output:
#        out = "VelocytoCounts/all.out"
#    params:
#        outf = "VelocytoCounts"
#    shell: """
#            export LC_ALL=en_US.utf-8
#            export LANG=en_US.utf-8
#            velocyto run --bcfile {input.bc} --outputfolder {params.outf} {input.bam} {input.gtf} > {output.out}
#    """
