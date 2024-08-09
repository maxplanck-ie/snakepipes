## function to get the name of the samplesheet and extend the name of the folder DESeq2 to DESeq2_[name]
def get_outdir(folder_name,sampleSheet):
    sample_name = os.path.splitext(os.path.basename(str(sampleSheet)))[0]

    return("{}_{}".format(folder_name, sample_name))

checkpoint split_sampleSheet:
    input:
        sampleSheet = sampleSheet
    output:
        splitSheets = os.path.join("splitSampleSheets",os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv")
    params:
        splitSheetPfx = os.path.join("splitSampleSheets",os.path.splitext(os.path.basename(str(sampleSheet)))[0])
    run:
        if isMultipleComparison:
            cf.splitSampleSheet(input.sampleSheet,params.splitSheetPfx)


## DESeq2 (on featureCounts)
rule DESeq2:
    input:
        counts_table = lambda wildcards : "featureCounts/counts_allelic.tsv" if 'allelic-mapping' in mode or "allelic-counting" in mode else "featureCounts/counts.tsv",
        sampleSheet = lambda wildcards: checkpoints.split_sampleSheet.get(compGroup=wildcards.compGroup).output,
        symbol_file = "Annotation/genes.filtered.symbol" #get_symbol_file
    output:
         "{}/DESeq2.session_info.txt".format(get_outdir("DESeq2",os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv"))
    benchmark:
        "{}/.benchmark/DESeq2.featureCounts.benchmark".format(get_outdir("DESeq2",os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv"))
    params:
        script=os.path.join(maindir, "shared", "rscripts", "DESeq2.R"),
        outdir = lambda wildcards,input: get_outdir("DESeq2",input.sampleSheet),
        sampleSheet = lambda wildcards,input: os.path.join(outdir,str(input.sampleSheet)),
        fdr = fdr,
        importfunc = os.path.join(maindir, "shared", "rscripts", "DE_functions.R"),
        allele_info = lambda wildcards : 'TRUE' if 'allelic-mapping' in mode or "allelic-counting" in mode else 'FALSE',
        tx2gene_file = 'NA',
        rmdTemplate = os.path.join(maindir, "shared", "rscripts", "DESeq2Report.Rmd"),
        formula = config["formula"],
        counts_table = lambda wildcards,input: os.path.join(outdir,input.counts_table),
        symbol_file = lambda wildcards,input: os.path.join(outdir,input.symbol_file)
    conda: CONDA_RNASEQ_ENV
    script: "{params.script}"

## DESeq2 (on Salmon)
rule DESeq2_Salmon_basic:
    input:
        counts_table = "Salmon/counts.transcripts.tsv",
        sampleSheet = lambda wildcards: checkpoints.split_sampleSheet.get(compGroup=wildcards.compGroup).output,
        tx2gene_file = "Annotation/genes.filtered.t2g",
        symbol_file = "Annotation/genes.filtered.symbol" #get_symbol_file
    output:
        "{}/DESeq2.session_info.txt".format(get_outdir("DESeq2_Salmon",os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv"))
    benchmark:
        "{}/.benchmark/DESeq2.Salmon.benchmark".format(get_outdir("DESeq2_Salmon",os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv"))
    params:
        script=os.path.join(maindir, "shared", "rscripts", "DESeq2.R"),
        outdir = lambda wildcards,input: get_outdir("DESeq2_Salmon",input.sampleSheet),
        sampleSheet = lambda wildcards,input: os.path.join(outdir,str(input.sampleSheet)),
        fdr = fdr,
        importfunc = os.path.join(maindir, "shared", "rscripts", "DE_functions.R"),
        allele_info = 'FALSE',
        tx2gene_file = os.path.join(outdir,"Annotation/genes.filtered.t2g"),
        rmdTemplate = os.path.join(maindir, "shared", "rscripts", "DESeq2Report.Rmd"),
        formula = config["formula"],
        counts_table = lambda wildcards,input: os.path.join(outdir,input.counts_table),
        symbol_file = lambda wildcards,input: os.path.join(outdir,input.symbol_file)
    conda: CONDA_RNASEQ_ENV
    script: "{params.script}"
