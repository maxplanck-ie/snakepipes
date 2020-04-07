import os

rule link_bam:
    input:
        indir + "/{sample}" + bamExt
    output:
        aligner + "/{sample}.unsorted.bam" if pipeline=="noncoding-rna-seq" else aligner + "/{sample}.bam"
    params:
        input_bai = indir + "/{sample}" + bamExt + ".bai",
        output_bai = aligner + "/{sample}.unsorted.bam.bai" if pipeline=="noncoding-rna-seq" else aligner + "/{sample}.bam.bai"
    run:
        if os.path.exists(params.input_bai) and not os.path.exists(os.path.join(outdir,params.output_bai)):
            os.symlink(params.input_bai,os.path.join(outdir,params.output_bai))
        if not os.path.exists(os.path.join(outdir,output[0])):
            os.symlink(os.path.join(outdir,input[0]),os.path.join(outdir,output[0]))

if not pipeline=="noncoding-rna-seq":
    rule samtools_index_external:
        input:
            aligner + "/{sample}.bam"
        output:
            aligner + "/{sample}.bam.bai"
        conda: CONDA_SHARED_ENV
        shell: "[[ ! -f {output} ]] samtools index {input}"

    if not pipeline=="WGBS" or pipeline=="WGBS" and skipBamQC:
        rule link_bam_bai_external:
            input:
                bam = aligner + "/{sample}.bam",
                bai = aligner + "/{sample}.bam.bai"
            output:
                bam_out = "filtered_bam/{sample}.filtered.bam",
                bai_out = "filtered_bam/{sample}.filtered.bam.bai",
            run:
                if not os.path.exists(os.path.join(outdir,output.bam_out)):
                    os.symlink(os.path.join(outdir,input.bam),os.path.join(outdir,output.bam_out))
                    os.symlink(os.path.join(outdir,input.bai),os.path.join(outdir,output.bai_out))


    rule sambamba_flagstat:
       input:
           aligner + "/{sample}.bam"
       output:
           "Sambamba/{sample}.markdup.txt"
       log: "Sambamba/logs/{sample}.flagstat.log"
       conda: CONDA_SAMBAMBA_ENV
       shell: """
           sambamba flagstat -p {input} > {output} 2> {log}
           """
