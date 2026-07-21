## function to get the name of the samplesheet and extend the name of the folder DESeq2 to DESeq2_[name]
def get_outdir(folder_name,sampleSheet,lrt):
    sample_name = os.path.splitext(os.path.basename(str(sampleSheet)))[0]
    output_folder_name = "{}_{}".format(folder_name, sample_name)
    if lrt:
        output_folder_name="{}_{}_LRT".format(folder_name, sample_name)
    return(output_folder_name)

## DESeq2 (on featureCounts)
rule DESeq2:
    input:
        counts_table = lambda wildcards : "featureCounts/counts_allelic.tsv" if 'allelic-mapping' in mode or 'allelic-counting' in mode or 'allelic-whatshap' in mode else "featureCounts/counts.tsv",
        sampleSheet = sampleSheet,
        symbol_file = "Annotation/genes.filtered.symbol" #get_symbol_file
    output:
        "{}/DESeq2.session_info.txt".format(get_outdir("DESeq2",sampleSheet,LRT)) 
    benchmark:
        "{}/.benchmark/DESeq2.featureCounts.benchmark".format(get_outdir("DESeq2",sampleSheet,LRT))
    params:
        script = os.path.join(maindir, "shared", "rscripts", "DESeq2.R"),
        sampleSheet = lambda wildcards,input: input.sampleSheet,
        outdir = get_outdir("DESeq2",sampleSheet,LRT),
        fdr = fdr,
        importfunc = os.path.join(maindir, "shared", "rscripts", "DE_functions.R"),
        allele_info = lambda wildcards : 'TRUE' if 'allelic-mapping' in mode or 'allelic-counting' in mode or 'allelic-whatshap' in mode else 'FALSE',
        tx2gene_file = 'NA',
        rmdTemplate = os.path.join(maindir, "shared", "rscripts", "DESeq2Report.Rmd"),
        formula = config["formula"],
        counts_table = lambda wildcards,input: os.path.join(outdir,input.counts_table),
        symbol_file = lambda wildcards,input: os.path.join(outdir,input.symbol_file),
        lrt = LRT
    conda: CONDA_RNASEQ_ENV
    script: "{params.script}"

## DESeq2 (on Salmon)
rule DESeq2_Salmon_basic:
    input:
        counts_table = "Salmon/counts.transcripts.tsv",
        sampleSheet = sampleSheet,
        tx2gene_file = "Annotation/genes.filtered.t2g",
        symbol_file = "Annotation/genes.filtered.symbol" #get_symbol_file
    output:
        "{}/DESeq2.session_info.txt".format(get_outdir("DESeq2_Salmon",sampleSheet,LRT))
    benchmark:
        "{}/.benchmark/DESeq2.Salmon.benchmark".format(get_outdir("DESeq2_Salmon",sampleSheet,LRT))
    params:
        script=os.path.join(maindir, "shared", "rscripts", "DESeq2.R"),
        sampleSheet = lambda wildcards,input: input.sampleSheet,
        outdir = get_outdir("DESeq2_Salmon",sampleSheet,LRT),
        fdr = fdr,
        importfunc = os.path.join(maindir, "shared", "rscripts", "DE_functions.R"),
        allele_info = 'FALSE',
        tx2gene_file = os.path.join(outdir,"Annotation/genes.filtered.t2g"),
        rmdTemplate = os.path.join(maindir, "shared", "rscripts", "DESeq2Report.Rmd"),
        formula = config["formula"],
        counts_table = lambda wildcards,input: os.path.join(outdir,input.counts_table),
        symbol_file = lambda wildcards,input: os.path.join(outdir,input.symbol_file),
        lrt = LRT
    conda: CONDA_RNASEQ_ENV
    script: "{params.script}"


rule DESeq2_Salmon_allelic:
    input:
        counts_table = "SalmonAllelic/counts.transcripts.tsv",
        sampleSheet = sampleSheet,
        tx2gene_file = "Annotation/genes.filtered.t2g",
        symbol_file = "Annotation/genes.filtered.symbol" #get_symbol_file
    output:
        "{}/DESeq2.session_info.txt".format(get_outdir("DESeq2_SalmonAllelic",sampleSheet,LRT))
    benchmark:
        "{}/.benchmark/DESeq2.SalmonAllelic.benchmark".format(get_outdir("DESeq2_SalmonAllelic",sampleSheet,LRT))
    params:
        script=os.path.join(maindir, "shared", "rscripts", "DESeq2.R"),
        sampleSheet = lambda wildcards,input: input.sampleSheet,
        outdir = get_outdir("DESeq2_SalmonAllelic",sampleSheet,LRT),
        fdr = fdr,
        importfunc = os.path.join(maindir, "shared", "rscripts", "DE_functions.R"),
        allele_info = 'TRUE',
        tx2gene_file = os.path.join(outdir,"Annotation/genes.filtered.t2g"),
        rmdTemplate = os.path.join(maindir, "shared", "rscripts", "DESeq2Report.Rmd"),
        formula = config["formula"],
        counts_table = lambda wildcards,input: os.path.join(outdir,input.counts_table),
        symbol_file = lambda wildcards,input: os.path.join(outdir,input.symbol_file),
        lrt = LRT
    conda: CONDA_RNASEQ_ENV
    script: "{params.script}"
