##STARsolo
##remember that reads are swapped in internals.snakefile!!
###currently having CB and UB tags output in the bam requires --outSAMtype SortedByCoordinate !!
import numpy
import os

rule STARsolo:
    input:
        r1="originalFASTQ/{sample}"+reads[0]+".fastq.gz",
        r2="originalFASTQ/{sample}"+reads[1]+".fastq.gz",
        annot="Annotation/genes.filtered.gtf"
    output:
        bam = "STARsolo/{sample}.sorted.bam",
        raw_counts = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/matrix.mtx",
        filtered_counts = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/matrix.mtx",
        filtered_bc = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/barcodes.tsv"
    params:
        alignerOptions = str(alignerOptions or ''),
        gtf = outdir+"/Annotation/genes.filtered.gtf",
        index = star_index,
        prefix = "STARsolo/{sample}/{sample}.",
        samsort_memory = '2G',
        sample_dir = "STARsolo/{sample}",
        bclist = BCwhiteList,
        UMIstart = STARsoloCoords[0],
        UMIlen = STARsoloCoords[1],
        CBstart = STARsoloCoords[2],
        CBlen = STARsoloCoords[3],
        outdir = outdir
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
            --outSAMunmapped Within \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS nM CB UB \
            --sjdbGTFfile {params.gtf} \
            --genomeDir {params.index} \
            --readFilesIn  {input.r1} {input.r2} \
            --readFilesCommand gunzip -c \
            --outFileNamePrefix {params.prefix} \
	    --soloType CB_UMI_Simple \
            --soloFeatures Gene Velocyto \
	    --soloUMIstart {params.UMIstart} \
	    --soloUMIlen {params.UMIlen} \
	    --soloCBstart {params.CBstart} \
	    --soloCBlen {params.CBlen} \
	    --soloCBwhitelist {params.bclist} \
	    --soloBarcodeReadLength 0 \
	    --soloCBmatchWLtype Exact \
	    --soloStrand Forward\
	    --soloUMIdedup Exact

        ln -s {params.outdir}/{params.prefix}Aligned.sortedByCoord.out.bam {params.outdir}/{output.bam}
 
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
           sambamba view -F "not unmapped and [CB] !=null" -t {threads} -f bam {input.bamfile} > {output.bamfile};
           sambamba index -t {threads} {output.bamfile}
           """

rule gzip_STARsolo_for_seurat:
    input:
        raw_counts = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/matrix.mtx",
        filtered_counts = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/matrix.mtx"
    output:
        raw_counts_gz = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/matrix.mtx.gz",
        filtered_counts_gz = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/matrix.mtx.gz"
    params:
        raw_bc = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/barcodes.tsv",
        filtered_bc = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/barcodes.tsv",
        raw_features = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/features.tsv",
        filtered_features = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/features.tsv",
        raw_bc_gz = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/barcodes.tsv.gz",
        filtered_bc_gz = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/barcodes.tsv.gz",
        raw_features_gz = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/features.tsv.gz",
        filtered_features_gz = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/features.tsv.gz"
    shell: """
         gzip -c {params.raw_bc} > {params.raw_bc_gz};
         gzip -c {params.raw_features} > {params.raw_features_gz};
         gzip -c {params.filtered_bc} > {params.filtered_bc_gz};
         gzip -c {params.filtered_features} > {params.filtered_features_gz};
         gzip -c {input.raw_counts} > {output.raw_counts_gz};
         gzip -c {input.filtered_counts} > {output.filtered_counts_gz}
    """


rule STARsolo_raw_to_seurat:
    input:
        infiles = expand("STARsolo/{sample}/{sample}.Solo.out/Gene/raw/matrix.mtx.gz",sample=samples)
    output:
        seurat = "Seurat/STARsolo_raw/merged_samples.RDS"
    params:
        indirs = expand(outdir + "/STARsolo/{sample}/{sample}.Solo.out/Gene/raw",sample=samples),
        wdir = outdir + "/Seurat/STARsolo_raw",
        samples = samples
    log:
        out = "Seurat/STARsolo_raw/logs/seurat.out"
    conda: CONDA_seurat3_ENV
    script: "../rscripts/scRNAseq_Seurat3.R"

rule STARsolo_filtered_to_seurat:
    input:
        infiles = expand("STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/matrix.mtx.gz",sample=samples)
    output:
        seurat = "Seurat/STARsolo_filtered/merged_samples.RDS"
    params:
        indirs = expand(outdir +"/STARsolo/{sample}/{sample}.Solo.out/Gene/filtered",sample=samples),
        wdir = outdir +"/Seurat/STARsolo_filtered",
        samples = samples
    log:
        out = "Seurat/STARsolo_filtered/logs/seurat.out"
    conda: CONDA_seurat3_ENV
    script: "../rscripts/scRNAseq_Seurat3.R"

if not skipVelocyto:
    rule cellsort_bam:
        input:
            bam = "filtered_bam/{sample}.filtered.bam"
        output:
            bam = "filtered_bam/cellsorted_{sample}.filtered.bam"
        params:
            samsort_memory="10G"
        threads: 4
        conda: CONDA_scRNASEQ_ENV
        shell: """
                MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX)
                samtools sort -m {params.samsort_memory} -@ {threads} -T $MYTEMP/{wildcards.sample} -t CB -O bam -o {output.bam} {input.bam}
                rm -rf $MYTEMP
               """

    #the barcode whitelist is currently taken from STARsolo filtered output, this is required to reduce runtime!
    #velocyto doesn't accept our filtered gtf; will have to use the mask, after all
    #no metadata table is provided

    checkpoint velocyto:
        input:
            gtf = "Annotation/genes.filtered.gtf",
            bam = "filtered_bam/{sample}.filtered.bam",
            csbam="filtered_bam/cellsorted_{sample}.filtered.bam",
            bc = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/barcodes.tsv"
        output:
            outdir = directory("VelocytoCounts/{sample}"),
            outdum = "VelocytoCounts/{sample}.done.txt"
        params:
            tempdir = tempDir
        conda: CONDA_scRNASEQ_ENV
        shell: """
                export LC_ALL=en_US.utf-8
                export LANG=en_US.utf-8
                export TMPDIR={params.tempdir}
                MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
                velocyto run --bcfile {input.bc} --outputfolder {output.outdir} --dtype uint64 {input.bam} {input.gtf};
                touch {output.outdum};
                rm -rf $MYTEMP
        """

    rule combine_loom:
        input: expand("VelocytoCounts/{sample}",sample=samples)
        output: "VelocytoCounts_merged/merged.loom"
        conda: CONDA_loompy_ENV
        params:
            outfile = outdir+"/VelocytoCounts_merged/merged.loom",
            script = maindir+"/shared/tools/loompy_merge.py",
            input_fp = lambda wildcards,input: [ os.path.join(outdir,f) for f in input ]
        shell: """
            python {params.script} -outf {params.outfile} {params.input_fp}
              """

    #rule velocity_to_seurat:
    #    input:
    #        indirs = expand("VelocytoCounts/{sample}",sample=samples)
    #    output:
    #        seurat = "Seurat/Velocyto/merged_samples.RDS"
    #    params:
    #        wdir = outdir + "/Seurat/Velocyto",
    #        samples = samples
    #    log:
    #        out = "Seurat/Velocyto/logs/seurat.out"
    #    conda: CONDA_seurat3_ENV
    #    script: "../rscripts/scRNAseq_merge_loom.R"
    
