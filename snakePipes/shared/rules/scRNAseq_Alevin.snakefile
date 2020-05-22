### cut columns 1 and 3 from t2g tsv #########

rule cut_t2g:
        input:
            "Annotation/genes.filtered.t2g"
        output:
            "Annotation/genes.filtered.slim.t2g"
        threads: 1
        shell:"""
            cut -f1,3 {input[0]} > {output[0]}
            """

rule SalmonAlevin:
        input:
            R2 = "originalFASTQ/{sample}"+reads[0]+".fastq.gz",
            R1 = "originalFASTQ/{sample}"+reads[1]+".fastq.gz",
            bin = "Salmon/SalmonIndex/seq.bin"
        params:
            protocol = "--" + prepprotocol,
            whitelist = "--whitelist {}".format(BCwhiteList) if BCwhiteList else "",
            expectcells = "--expectcells {}".format(expectcells) if expectcells else "",
            tgMap = "Annotation/genes.filtered.slim.t2g",
            libtype = libType,
            outdir = "Alevin/{sample}"
        output:
            quantmat = "Alevin/{sample}/alevin/quants_mat.gz"
        threads: 40
        shell:"""
            salmon alevin -l {params.libtype} -1 {input.R1} -2 {input.R2} {params.protocol} -i Salmon/SalmonIndex -p {threads} -o {params.outdir} --tgMap {params.tgMap} --dumpFeatures --dumpMtx --numCellBootstraps 100
            """
rule AlevinQC:
        input:
            indir = expand("Alevin/{sample}/", sample=samples)
        output:
            outfiles = expand("QC/Alevin_{sample}.html", sample=samples)
        params:
            script=os.path.join(maindir, "shared", "rscripts", "scRNAseq_Alevinqc.R"),
            outdir =  "QC/",
            samid = expand("{sample}", sample=samples),
            outfile = expand("Alevin_{sample}.html", sample=samples)
        log:
            out = expand("Alevin/{sample}/QC.out", sample=samples),
        conda: CONDA_scRNASEQ_ENV
        script: """
            Rscript {params.script} {input.indir} {params.outdir} {params.samid} {params.outfile} >& {log.out}
            echo "{params.script} {input.indir} {params.outdir} {params.samid} {params.outfile}"
            """
