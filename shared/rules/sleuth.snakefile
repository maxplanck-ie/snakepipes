## sleuth (on Salmon)
rule sleuth_Salmon:
    input:
        quant_files = expand("Salmon/{sample}/abundance.h5", sample=samples),
        t2g = 'Annotation/genes.filtered.t2g',
        sample_info = sample_info
    output:
        "sleuth/so.rds"
    benchmark:
        "sleuth/.benchmark/sleuth.Salmon.benchmark"
    params:
        indir = "Salmon",
        outdir = "sleuth",
        fdr = 0.05,
    log: "sleuth/sleuth.log"
    conda: CONDA_RNASEQ_ENV
    shell:
#        "( cd {params.outdir} && export R_LIBS_USER="+R_libs_path+" && "
#        "cat "+os.path.join(workflow_tools,"sleuth.R")+" | "
#        ""+os.path.join(R_path,"R")+" --vanilla --args "
        "cd {param.oudir} && "
        "Rscipt "+os.path.join(maindir, "shared", "tools", "sleuth.R")
        "{input.sample_info} "
        "../{params.indir} "
        "../{params.outdir} "
        "{params.fdr} "
        "../{input.t2g} "
        "2>&1 | tee {log}"

## sleuth (on Salmon)
