## function to get the name of the samplesheet and extend the name of the folder DESeq2 to DESeq2_[name]
def get_outdir(folder_name,sampleSheet):
    sample_name = os.path.splitext(os.path.basename(str(sampleSheet)))[0]

    return("{}_{}".format(folder_name, sample_name))

## DESeq2 (on featureCounts)
rule DESeq2:
    input:
        counts_table = lambda wildcards : "featureCounts/counts_allelic.tsv" if 'allelic-mapping' in mode else "featureCounts/counts.tsv",
        sampleSheet = sampleSheet,
        symbol_file = "Annotation/genes.filtered.symbol" #get_symbol_file
    output:
        "{}/DESeq2.session_info.txt".format(get_outdir("DESeq2",sampleSheet))
    benchmark:
        "{}/.benchmark/DESeq2.featureCounts.benchmark".format(get_outdir("DESeq2",sampleSheet))
    params:
        script=os.path.join(maindir, "shared", "rscripts", "DESeq2.R"),
        outdir = get_outdir("DESeq2",sampleSheet),
        fdr = 0.05,
        importfunc = os.path.join(maindir, "shared", "rscripts", "DE_functions.R"),
        allele_info = lambda wildcards : 'TRUE' if 'allelic-mapping' in mode else 'FALSE',
        tx2gene_file = 'NA',
        rmdTemplate = os.path.join(maindir, "shared", "rscripts", "DESeq2Report.Rmd")
    log:
        out = "{}/logs/DESeq2.out".format(get_outdir("DESeq2",sampleSheet)),
        err = "{}/logs/DESeq2.err".format(get_outdir("DESeq2",sampleSheet))
    conda: CONDA_RNASEQ_ENV
    shell:
        "cd {params.outdir} && "
        "Rscript {params.script} "
        "{input.sampleSheet} " # 1
        "../{input.counts_table} " # 2
        "{params.fdr} " # 3
        "../{input.symbol_file} " # 4
        "{params.importfunc} " # 5
        "{params.allele_info} " # 6
        "{params.tx2gene_file} " # 7
        "{params.rmdTemplate} " # 8
        " > ../{log.out} 2> ../{log.err}"


## DESeq2 (on Salmon)
rule DESeq2_Salmon:
    input:
        counts_table = "Salmon/counts.transcripts.tsv" if not "allelic-mapping" in mode else "SalmonAllelic/counts.transcripts.tsv",
        sampleSheet = sampleSheet,
        tx2gene_file = "Annotation/genes.filtered.t2g",
        symbol_file = "Annotation/genes.filtered.symbol" #get_symbol_file
    output:
        "{}/DESeq2.session_info.txt".format(get_outdir("DESeq2_Salmon",sampleSheet))
    log:
        out = "{}/logs/DESeq2.out".format(get_outdir("DESeq2_Salmon",sampleSheet)),
        err = "{}/logs/DESeq2.err".format(get_outdir("DESeq2_Salmon",sampleSheet))
    benchmark:
        "{}/.benchmark/DESeq2.Salmon.benchmark".format(get_outdir("DESeq2_Salmon",sampleSheet))
    params:
        script=os.path.join(maindir, "shared", "rscripts", "DESeq2.R"),
        outdir = get_outdir("DESeq2_Salmon",sampleSheet),
        fdr = 0.05,
        importfunc = os.path.join(maindir, "shared", "rscripts", "DE_functions.R"),
        allele_info = 'FALSE',
        tx2gene_file = "Annotation/genes.filtered.t2g",
        rmdTemplate = os.path.join(maindir, "shared", "rscripts", "DESeq2Report.Rmd")
    conda: CONDA_RNASEQ_ENV
    shell:
        "cd {params.outdir} && "
        "Rscript {params.script} "
        "{input.sampleSheet} " # 1
        "../{input.counts_table} " # 2
        "{params.fdr} " # 3
        "../{input.symbol_file} " # 4
        "{params.importfunc} " # 5
        "{params.allele_info} " # 6
        "../{input.tx2gene_file} " # 7
        "{params.rmdTemplate} " # 8
        " > ../{log.out} 2> ../{log.err}"
