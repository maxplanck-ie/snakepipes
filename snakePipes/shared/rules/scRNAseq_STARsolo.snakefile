##STARsolo
##remember that reads are swapped in internals.snakefile!!
###currently having CB and UB tags output in the bam requires --outSAMtype SortedByCoordinate !!
rule STARsolo:
    input:
        r1="originalFASTQ/{sample}"+reads[0]+".fastq.gz",
        r2="originalFASTQ/{sample}"+reads[1]+".fastq.gz",
        annot="Annotation/genes.filtered.gtf"
    output:
        bam = "STARsolo/{sample}.sorted.bam"
    params:
        alignerOptions = str(alignerOptions or ''),
        gtf = outdir+"/Annotation/genes.filtered.gtf",
        index = star_index,
        prefix = "STARsolo/{sample}/{sample}.",
        samsort_memory = '2G',
        sample_dir = "STARsolo/{sample}",
        bam = "{sample}/{sample}.Aligned.sortedByCoord.out.bam"
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
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS nM CB UB \
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
	    --soloUMIdedup Exact ;

        ln -s {params.bam} {output.bam};
 
        rm -rf $MYTEMP
         """


rule filter_bam:
    input:
        bamfile = aligner+"/{sample}.sorted.bam",
        bami = aligner+"/{sample}.sorted.bam.bai"
    output:
        bamfile = "filtered_bam/{sample}.filtered.bam",
        bami = "filtered_bam/{sample}.filtered.bam.bai"
    threads: 8
    conda: CONDA_SAMBAMBA_ENV
    shell: """
           MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
           sambamba view -F "not unmapped and [CB] !=null" -t {threads} -f bam {input.bamfile} > {output.bamfile};
           samtools index {output.bamfile};
           rm -rf $MYTEMP
           """

rule cellsort_bam:
    input:
        bam = "filtered_bam/{sample}.filtered.bam"
    output:
        bam = "filtered_bam/cellsorted_{sample}.bam"
#        bami = "filtered_bam/cellsorted_{sample}.bam.bai"
    params:
        samsort_memory="3G"
    threads: 4
    conda: CONDA_scRNASEQ_ENV
    shell: """
            MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
            samtools sort -m {params.samsort_memory} -@ {threads} -T $MYTEMP/{wildcards.sample} -t CB -O bam -o {output.bam} {input.bam};
            rm -rf $MYTEMP
           """

#the barcode whitelist is currently passed in although it's not tested if it's actually necessery as it was already provided to STARsolo
#gtf mask is not used as a filtered gtf is passed in
#no metadata table is provided

rule velocyto:
    input:
        bc = "/data/processing/bioinfo-core/celseq_barcodes.384.1col.txt",
        gtf = "Annotation/genes.filtered.gtf",
        bam = "filtered_bam/cellsorted_{sample}.bam"
    output:
        out = "VelocytoCounts/{sample}.loom"
    params:
        outf = "VelocytoCounts"
    shell: """
            export LC_ALL=en_US.utf-8
            export LANG=en_US.utf-8
            velocyto run --bcfile {input.bc} --outputfolder {params.outf} {input.bam} {input.gtf} 
    """
