## sleuth (on Salmon)
rule sleuth_prep:
    input:
        genes_gtf
    output:
        # DEseq2 path
        t2g
    benchmark:
        "sleuth/.benchmark/sleuth_prep.Salmon.benchmark"
    log: "sleuth/sleuth_prep.log"
    shell:
        os.path.join(workflow_tools, 'extractTx2GeneMapping.sh') + "{input} {output}"
