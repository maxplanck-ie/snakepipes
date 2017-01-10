## If a symbol file is available
def get_symbol_file(wildcards):
    symbol_file = os.path.join(maindir, "shared", "organisms", genome+".symbol")
    if os.path.isfile(symbol_file):
        return(symbol_file)
    else:
        return("")


## DESeq2 (on featureCounts)
rule DESeq2:
    input:
        counts_table = "featureCounts/counts.tsv",
        sample_info = sample_info
    output:
        "DESeq2/DESeq2.session_info.txt"
    benchmark:
        "DESeq2/.benchmark/DESeq2.featureCounts.benchmark"
    params:
        outdir = "DESeq2",
        fdr = 0.05,
        symbol_file = get_symbol_file
    log: "DESeq2/DESeq2.log"
    shell:
        "( cd {params.outdir} && export R_LIBS_USER="+R_libs_path+" && "
        "cat "+os.path.join(workflow_tools,"DESeq2.R")+" | "
        ""+os.path.join(R_path,"R")+" --vanilla --args "
        "{input.sample_info} "
        "../{input.counts_table} "
        "{params.fdr} "
        "{params.symbol_file} "
        ") 2>&1 | tee {log}"


## DESeq2 (on Salmon)
rule DESeq2_Salmon:
    input:
        counts_table = "Salmon/counts.genes.tsv",
        sample_info = sample_info
    output:
        "DESeq2_Salmon/DESeq2.session_info.txt"
    benchmark:
        "DESeq2_Salmon/.benchmark/DESeq2.Salmon.benchmark"
    params:
        outdir = "DESeq2_Salmon",
        fdr = 0.05,
        symbol_file = get_symbol_file
    log: "DESeq2_Salmon/DESeq2.log"
    shell:
        "( cd {params.outdir} && export R_LIBS_USER="+R_libs_path+" && "
        "cat "+os.path.join(workflow_tools,"DESeq2.R")+" | "
        ""+os.path.join(R_path,"R")+" --vanilla --args "
        "{input.sample_info} "
        "../{input.counts_table} "
        "{params.fdr} "
        "{params.symbol_file} "
        ") 2>&1 | tee {log}"
