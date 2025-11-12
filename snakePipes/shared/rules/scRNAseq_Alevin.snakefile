### cut columns 1 and 3 from t2g tsv #########
import gzip
import os

def get_flank_length(file,read_length_frx):
    if not os.path.exists(file):
        return (0)
    with gzip.open(file,"r") as rf:
        head = [next(rf) for x in range(4)]
    read_length = len(head[1])
    flank_length = int( (1 - read_length_frx) * read_length )
    return(flank_length)

rule cut_t2g:
        input:
            genes_t2g
        output:
            "Annotation/genes.slim.t2g"
        threads: 1
        shell:"""
            cut -f1,2 {input[0]} > {output[0]}
            """

rule SalmonAlevin:
        input:
            R2 = "originalFASTQ/{sample}"+reads[0]+".fastq.gz",
            R1 = "originalFASTQ/{sample}"+reads[1]+".fastq.gz",
            tgMap = "Annotation/genes.slim.t2g"
        params:
            index = salmon_index,
            protocol = "--" + prepProtocol,
            whitelist = "--whitelist {}".format(BCwhiteList) if BCwhiteList else "",
            expectcells = "--expectcells {}".format(expectCells) if expectCells else "",
            libtype = alevinLibraryType,
            outdir = "Alevin/{sample}"
        output:
            quantmat = "Alevin/{sample}/alevin/quants_mat.gz",
        #Use RNAseq env because Salmon already installed there (no need for duplication).
        conda: CONDA_SALMON_ENV
        threads: 8
        shell:"""
            salmon alevin -l {params.libtype} -1 {input.R1} -2 {input.R2} {params.protocol} -i {params.index} -p {threads} -o {params.outdir} --tgMap {input.tgMap} --dumpFeatures --dumpMtx --numCellBootstraps 100
            """

rule AlevinQC:
        input:
            indum = "Alevin/{sample}/alevin/quants_mat.gz"
        output:
            outfiles = "multiQC/Alevin_{sample}.html"
        params:
            indir = "Alevin/{sample}/",
            outdir =  "multiQC/",
            samid = "{sample}",
            outfile = "Alevin_{sample}.html"
        conda: CONDA_alevinqc_ENV
        script: "../rscripts/scRNAseq_Alevinqc.R"


rule AlevinForVelocity:
        input:
            R2 = "originalFASTQ/{sample}"+reads[0]+".fastq.gz",
            R1 = "originalFASTQ/{sample}"+reads[1]+".fastq.gz"
        params:
            velo_index = salmon_velocity_index,
            tgMap = t2g_velocity,
            protocol = "--" + prepProtocol,
            whitelist = "--whitelist {}".format(BCwhiteList) if BCwhiteList else "",
            expectcells = "--expectcells {}".format(expectCells) if expectCells else "",
            libtype = alevinLibraryType,
            outdir = "AlevinForVelocity/{sample}"
        output:
            quantmat = "AlevinForVelocity/{sample}/alevin/quants_mat.gz"
        #Use RNAseq env because Salmon already installed there (no need for duplication).
        conda: CONDA_SALMON_ENV
        threads: 8
        shell:"""
            salmon alevin -l {params.libtype} -1 {input.R1} -2 {input.R2} {params.protocol} -i {params.velo_index} -p {threads} -o {params.outdir} --tgMap {params.tgMap} --dumpFeatures --dumpMtx --numCellBootstraps 100
            """

rule velo_to_sce:
    input:
        quantmat = expand("AlevinForVelocity/{sample}/alevin/quants_mat.gz",sample=samples),
        t2g = t2g_velocity,
        g2s = "Annotation/genes.filtered.symbol"
    output:
        merged = "SingleCellExperiment/AlevinForVelocity/merged_samples.RDS"
    params:
        wdir = os.path.join(outdir,"SingleCellExperiment/AlevinForVelocity"),
        alevindir = os.path.join(outdir,"AlevinForVelocity"),
        samplenames = samples,
        t2g = lambda wildcards,input: os.path.join(outdir, input.t2g),
        g2s = lambda wildcards,input: os.path.join(outdir, input.g2s),
        outfile = lambda wildcards,output: os.path.join(outdir, output.merged)
    conda: CONDA_eisaR_ENV
    script: "../rscripts/scRNAseq_splitAlevinVelocityMatrices.R"
