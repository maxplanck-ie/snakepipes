## sleuth (on Salmon)

sample_name = os.path.splitext(os.path.basename(sampleSheet))[0]

rule sleuth_Salmon:
    input:
        quant_files = expand("Salmon/{sample}/abundance.h5", sample=samples),
        t2g = "Annotation/genes.filtered.t2g",
        sampleSheet = sampleSheet
    output:
        "sleuth_Salmon_{}/so.rds".format(sample_name)
    benchmark:
        "sleuth/.benchmark/sleuth.Salmon.benchmark"
    params:
        script=os.path.join(maindir, "shared", "rscripts", "sleuth.R"),
        indir = os.path.join(outdir,"Salmon"),
        outdir = os.path.join(outdir,"sleuth_Salmon_{}".format(sample_name)),
        fdr = 0.05,
    log:
        out = "sleuth_Salmon_{}/logs/sleuth.out".format(sample_name),
        err = "sleuth_Salmon_{}/logs/sleuth.err".format(sample_name)
    conda: CONDA_RNASEQ_ENV
    shell:
        "Rscript {params.script} "
        "{input.sampleSheet} "
        "{params.indir} "
        "{params.outdir} "
        "{params.fdr} " + os.path.join(outdir,"{input.t2g}") +
        " >{log.out} 2>{log.err}"
