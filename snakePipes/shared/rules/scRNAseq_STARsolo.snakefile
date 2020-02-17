##STARsolo
##remember that reads are swapped in internals.snakefile!!
###currently having CB and UB tags output in the bam requires --outSAMtype SortedByCoordinate !!
import numpy
#import loompy
import os

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
        bam = "{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bclist = BCwhiteList,
        UMIstart = STARsoloCoords[0],
        UMIlen = STARsoloCoords[1],
        CBstart = STARsoloCoords[2],
        CBlen = STARsoloCoords[3]
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
           sambamba view -F "not unmapped and [CB] !=null" -t {threads} -f bam {input.bamfile} > {output.bamfile};
           sambamba index -t {threads} {output.bamfile}
           """

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
            MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
            samtools sort -m {params.samsort_memory} -@ {threads} -T $MYTEMP/{wildcards.sample} -t CB -O bam -o {output.bam} {input.bam};
            rm -rf $MYTEMP
           """

#the barcode whitelist is currently passed in although it's not tested if it's actually necessery as it was already provided to STARsolo
#velocyto doesn't accept our filtered gtf; will have to use the mask, after all
#no metadata table is provided

checkpoint velocyto:
    input:
        gtf = genes_gtf,
        bam = "filtered_bam/{sample}.filtered.bam",
        csbam="filtered_bam/cellsorted_{sample}.filtered.bam"
    output:
        outdir = directory("VelocytoCounts/{sample}"),
        outdum = "VelocytoCounts/{sample}.done.txt"
    params:
        bc = BCwhiteList_gz,
        tempdir = tempDir
    conda: CONDA_scRNASEQ_ENV
    shell: """
            export LC_ALL=en_US.utf-8
            export LANG=en_US.utf-8
            export TMPDIR={params.tempdir}
            MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
            velocyto run --bcfile {params.bc} --outputfolder {output.outdir} --dtype uint64 {input.bam} {input.gtf};
            touch {output.outdum};
            rm -rf $MYTEMP
    """

#rule combine_loom:
#    input: expand("VelocytoCounts/{sample}",sample=samples)
#    output: "VelocytoCounts_merged/merged.txt"
#    run: 
#        filelist=[]
#        for p in input:
#            z=os.listdir(p)
#            f=list(filter(lambda x: '.loom' in x,z))
#            ifi=os.path.join(outdir,p,f[0])
#            filelist.append(ifi)
#        print(filelist)
#        outf=outdir+"/VelocytoCounts_merged/merged.loom"
#        print(outf)
#        loompy.combine(files=filelist,output_file=outf, key="Accession")
