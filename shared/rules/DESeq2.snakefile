## If a symbol file is available
#def get_symbol_file(wildcards):
#    symbol_file = os.path.join(maindir, "shared", "organisms", genome+".symbol")
#    if os.path.isfile(symbol_file):
#        return(symbol_file)
#    else:
#        return("")

## function to get the name of the samplesheet and extend the name of the folder DESeq2 to DESeq2_[name]
def get_outdir(folder_name):
    sample_name = os.path.splitext(os.path.basename(str(sample_info)))[0]
    return("{}_{}".format(folder_name, sample_name))

## DESeq2 (on featureCounts)
rule DESeq2:
    input:
        counts_table = lambda wildcards : "featureCounts/counts_allelic.tsv" if 'allelic-mapping' in mode else "featureCounts/counts.tsv",
        sample_info = sample_info,
        symbol_file = "Annotation/genes.filtered.symbol" #get_symbol_file
    output:
        "{}/DESeq2.session_info.txt".format(get_outdir("DESeq2"))
    benchmark:
        "{}/.benchmark/DESeq2.featureCounts.benchmark".format(get_outdir("DESeq2"))
    params:
        outdir = get_outdir("DESeq2"),
        fdr = 0.05,
        importfunc = os.path.join(workflow_tools,"snakediff","R", "DE_functions.R"),
        allele_info = lambda wildcards : 'TRUE' if 'allelic-mapping' in mode else 'FALSE',
        tx2gene_file = 'NA'
    log: "{}/DESeq2.log".format(get_outdir("DESeq2"))
    shell:
        "( cd {params.outdir} && export R_LIBS_USER="+R_libs_path+" && "
        "cat "+os.path.join(workflow_tools,"DESeq2.R")+" | "
        ""+os.path.join(R_path,"R")+" --vanilla --slave --args "
        "{input.sample_info} " # 1
        "../{input.counts_table} " # 2
        "{params.fdr} " # 3
        "../{input.symbol_file} " # 4
        "{params.importfunc} " # 5
        "{params.allele_info} " # 6
        "{params.tx2gene_file} " # 7
        " ) 2>&1 | tee {log}"


## DESeq2 (on Salmon)
rule DESeq2_Salmon:
    input:
        counts_table = "Salmon/counts.tsv",
        sample_info = sample_info,
        tx2gene_file = "Annotation/genes.filtered.t2g",
        symbol_file = "Annotation/genes.filtered.symbol" #get_symbol_file
    output:
        "{}/DESeq2.session_info.txt".format(get_outdir("DESeq2_Salmon"))
    benchmark:
        "{}/.benchmark/DESeq2.Salmon.benchmark".format(get_outdir("DESeq2_Salmon"))
    params:
        outdir = get_outdir("DESeq2_Salmon"),
        fdr = 0.05,
        importfunc = os.path.join(workflow_tools,"snakediff", "R" ,"DE_functions.R"),
        allele_info = 'FALSE',
        tx2gene_file = "Annotation/genes.filtered.t2g"
    log: "{}/DESeq2.log".format(get_outdir("DESeq2_Salmon"))
    shell:
        "( cd {params.outdir} && export R_LIBS_USER="+R_libs_path+" && "
        "cat "+os.path.join(workflow_tools,"DESeq2.R")+" | "
        ""+os.path.join(R_path,"R")+" --vanilla --slave --args "
        "{input.sample_info} " # 1
        "../{input.counts_table} " # 2
        "{params.fdr} " # 3
        "../{input.symbol_file} " # 4
        "{params.importfunc} " # 5
        "{params.allele_info} " # 6
        "../{input.tx2gene_file} " # 7
        ") 2>&1 | tee {log}"
