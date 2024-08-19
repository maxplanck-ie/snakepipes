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
        filtered_bc = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/barcodes.tsv",
        raw_features = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/features.tsv",
        filtered_features = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/features.tsv",
        summary = "STARsolo/{sample}/{sample}.Solo.out/Gene/Summary.csv"
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
        outdir = outdir,
        tempDir = tempDir
    benchmark:
        aligner+"/.benchmark/STARsolo.{sample}.benchmark"
    threads: lambda wildcards: 20 if 20<max_thread else max_thread  # 3.2G per core
    conda: CONDA_scRNASEQ_ENV
    shell: """
        TMPDIR={params.tempDir}
        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
        ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
        STAR --runThreadN {threads} \
            {params.alignerOptions} \
            --sjdbOverhang 100 \
            --outSAMunmapped Within \
            --outSAMtype BAM SortedByCoordinate \
            --outBAMsortingBinsN 20 \
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

        ln -s ../{params.prefix}Aligned.sortedByCoord.out.bam {output.bam}

        rm -rf $MYTEMP
         """

rule STARsolo_report:
    input:  expand("STARsolo/{sample}/{sample}.Solo.out/Gene/Summary.csv",sample=samples)
    output:
        report = "STARsolo/Report.tsv"
    params:
        wdir = outdir + "/STARsolo",
        input = lambda wildcards,input: [ os.path.join(outdir,x) for x in input ],
        samples = samples
    conda: CONDA_seurat_ENV
    script: "../rscripts/scRNAseq_report.R"


rule filter_bam:
    input:
        bamfile = aligner+"/{sample}.sorted.bam",
        bami = aligner+"/{sample}.sorted.bam.bai"
    output:
        bamfile = "filtered_bam/{sample}.filtered.bam",
        bami = "filtered_bam/{sample}.filtered.bam.bai"
    threads: lambda wildcards: 8 if 8<max_thread else max_thread
    conda: CONDA_SAMBAMBA_ENV
    shell: """
           sambamba view -F "not unmapped and [CB] !=null" -t {threads} -f bam {input.bamfile} > {output.bamfile};
           sambamba index -t {threads} {output.bamfile}
           """

##remove this rule as soon as STARsolo output has been fixed by Alex Dobin
rule STARsolo_features_to_V3:
    input:
        raw_features = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/features.tsv",
        filtered_features = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/features.tsv"
    output:
        raw_features = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/features.v3.tsv",
        filtered_features = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/features.v3.tsv"
    shell: """
        awk '{{print $1, $2, "."}}' {input.raw_features} | tr " " "\t" > {output.raw_features};
        awk '{{print $1, $2, "."}}' {input.filtered_features} | tr " " "\t" > {output.filtered_features}
    """


rule gzip_STARsolo_for_seurat:
    input:
        raw_counts = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/matrix.mtx",
        filtered_counts = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/matrix.mtx",
        raw_features = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/features.v3.tsv",
        filtered_features = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/features.v3.tsv"
    output:
        raw_counts_gz = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/matrix.mtx.gz",
        filtered_counts_gz = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/matrix.mtx.gz"
    params:
        raw_bc = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/barcodes.tsv",
        filtered_bc = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/barcodes.tsv",
        raw_bc_gz = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/barcodes.tsv.gz",
        filtered_bc_gz = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/barcodes.tsv.gz",
        raw_features_gz = "STARsolo/{sample}/{sample}.Solo.out/Gene/raw/features.tsv.gz",
        filtered_features_gz = "STARsolo/{sample}/{sample}.Solo.out/Gene/filtered/features.tsv.gz"
    shell: """
         gzip -c {params.raw_bc} > {params.raw_bc_gz};
         gzip -c {input.raw_features} > {params.raw_features_gz};
         gzip -c {params.filtered_bc} > {params.filtered_bc_gz};
         gzip -c {input.filtered_features} > {params.filtered_features_gz};
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
    conda: CONDA_seurat_ENV
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
    conda: CONDA_seurat_ENV
    script: "../rscripts/scRNAseq_Seurat3.R"




rule remove_empty_drops:
    input:
        infiles = expand("STARsolo/{sample}/{sample}.Solo.out/Gene/raw/matrix.mtx.gz",sample=samples)
    output:
        seurat = "Seurat/STARsolo_raw_RmEmptyDrops/merged_samples.RDS"
    params:
        indirs = expand(outdir + "/STARsolo/{sample}/{sample}.Solo.out/Gene/raw",sample=samples),
        wdir = outdir + "/Seurat/STARsolo_raw_RmEmptyDrops",
        samples = samples
    conda: CONDA_seurat_ENV
    script: "../rscripts/scRNAseq_EmptyDrops.R"


if not skipVelocyto:
    rule cellsort_bam:
        input:
            bam = "filtered_bam/{sample}.filtered.bam"
        output:
            bam = "filtered_bam/cellsorted_{sample}.filtered.bam"
        params:
            samsort_memory="10G",
            tempDir = tempDir
        threads: lambda wildcards: 4 if 4<max_thread else max_thread
        conda: CONDA_scRNASEQ_ENV
        shell: """
                TMPDIR={params.tempDir}
                MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX)
                samtools sort -m {params.samsort_memory} -@ {threads} -T $MYTEMP/{wildcards.sample} -t CB -O bam -o {output.bam} {input.bam}
                rm -rf $MYTEMP
               """

    #the barcode whitelist is currently taken from STARsolo filtered output, this is required to reduce runtime!
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
                velocyto run --bcfile {input.bc} --outputfolder {output.outdir} --dtype uint64 {input.bam} {input.gtf} ;
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
    #    conda: CONDA_seurat3_ENV
    #    script: "../rscripts/scRNAseq_merge_loom.R"
