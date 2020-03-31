import os
rule link_bam:
    input:
        indir + "/{sample}" + bamExt
    output:
        aligner + "/{sample}.unsorted.bam" if pipeline=="noncoding-rna-seq" else aligner + "/{sample}.bam"
    run:
        if not os.path.exists(os.path.join(outdir,output)):
            os.symlink(os.path.join(outdir,input),os.path.join(outdir,output))

if not pipeline=="noncoding-rna-seq":
    rule samtools_index_external:
        input:
            aligner + "/{sample}.bam"
        output:
            aligner + "/{sample}.bam.bai"
        conda: CONDA_SHARED_ENV
        shell: "samtools index {input}"


    rule link_bam_bai_external:
        input:
            bam = aligner + "/{sample}.bam",
            bai = aligner + "/{sample}.bam.bai"
        output:
            bam_out = "filtered_bam/{sample}.filtered.bam",
            bai_out = "filtered_bam/{sample}.filtered.bam.bai",
        run:
            if not os.path.exists(os.path.join(outdir,bam_out)):
                os.symlink(os.path.join(outdir,bam),os.path.join(outdir,bam_out))
                os.symlink(os.path.join(outdir,bai),os.path.join(outdir,bai_out))


    rule sambamba_flagstat:
       input:
           aligner + "/{sample}.bam"
       output:
           "Sambamba/{sample}.markdup.txt"
       conda: CONDA_SAMBAMBA_ENV
       shell: """
           sambamba flagstat -p {input} > {output}
           """
