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
            quantmat = "Alevin/{sample}/alevin/quants_mat.gz",
        log: 
            out =  "Alevin/logs/alevin.{sample}.out",
            err = "Alevin/logs/alevin.{sample}.err"
        #Use RNAseq env because Salmon already installed there (no need for duplication).
        conda: CONDA_RNASEQ_ENV
        threads: 40
        shell:"""
            salmon alevin -l {params.libtype} -1 {input.R1} -2 {input.R2} {params.protocol} -i Salmon/SalmonIndex -p {threads} -o {params.outdir} --tgMap {params.tgMap} --dumpFeatures --dumpMtx --numCellBootstraps 100 > {log.out} 2> {log.err}
            """
rule AlevinQC:
        input:
            indum = "Alevin/{sample}/alevin/quants_mat.gz"
        output:
            outfiles = "multiQC/Alevin_{sample}.html"
        params:
            indir = "Alevin/{sample}/",
            script = os.path.join(maindir, "shared", "rscripts", "scRNAseq_Alevinqc.R"),
            outdir =  "multiQC/",
            samid = "{sample}",
            outfile = "Alevin_{sample}.html"
        conda: CONDA_alevinqc_ENV
        shell: """
            Rscript {params.script} {params.indir} {params.outdir} {params.samid} {params.outfile}
            """
