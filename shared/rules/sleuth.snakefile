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
        script=os.path.join(maindir, "shared", "tools", "sleuth.R"),
        indir = "Salmon",
        outdir = "sleuth",
        fdr = 0.05,
    log: "sleuth.log"
    conda: CONDA_RNASEQ_ENV
    shell:
        "cd {params.outdir} && "
        "Rscript {params.script} "
        "{input.sample_info} "
        "../{params.indir} "
        "../{params.outdir} "
        "{params.fdr} "
        "../{input.t2g} "
        "2>&1 | tee {log}"

## sleuth (on Salmon)
