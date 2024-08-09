import os

if pipeline=="rnaseq" and "allelic-counting" in mode:
    rule link_bam:
        input:
            indir + "/{sample}.{suffix}" + bamExt
        output:
            "allelic_bams/{sample}.{suffix}" + bamExt
        params:
            input_bai = indir + "/{sample}.{suffix}" + bamExt + ".bai",
            output_bai = "allelic_bams/{sample}.{suffix}" + bamExt + ".bai"
        run:
            if os.path.exists(params.input_bai) and not os.path.exists(os.path.join(outdir,params.output_bai)):
                os.symlink(params.input_bai,os.path.join(outdir,params.output_bai))
            if not os.path.exists(os.path.join(outdir,output[0])):
                os.symlink(os.path.join(outdir,input[0]),os.path.join(outdir,output[0]))

    rule samtools_index_external:
        input:
            "allelic_bams/{sample}.{suffix}" + bamExt
        output:
            "allelic_bams/{sample}.{suffix}" + bamExt + ".bai"
        conda: CONDA_SHARED_ENV
        shell: "if [[ ! -f {output[0]} ]]; then samtools index {input[0]}; fi"


else:
    rule link_bam:
        input:
            indir + "/{sample}" + bamExt
        output:
            aligner + "/{sample}.unsorted.bam" if pipeline=="ncRNAseq" else aligner + "/{sample}.bam"
        params:
            input_bai = indir + "/{sample}" + bamExt + ".bai",
            output_bai = aligner + "/{sample}.unsorted.bam.bai" if pipeline=="ncRNAseq" else aligner + "/{sample}.bam.bai"
        run:
            if os.path.exists(params.input_bai) and not os.path.exists(os.path.join(outdir,params.output_bai)):
                os.symlink(params.input_bai,os.path.join(outdir,params.output_bai))
            if not os.path.exists(os.path.join(outdir,output[0])):
                os.symlink(os.path.join(outdir,input[0]),os.path.join(outdir,output[0]))

    if not pipeline=="ncRNAseq":
        rule samtools_index_external:
            input:
                aligner + "/{sample}.bam"
            output:
                aligner + "/{sample}.bam.bai"
            conda: CONDA_SHARED_ENV
            shell: "if [[ ! -f {output[0]} ]]; then samtools index {input[0]}; fi"

        if not pipeline=="WGBS" or pipeline=="WGBS" and skipBamQC:
            rule link_bam_bai_external:
                input:
                    bam = aligner + "/{sample}.bam",
                    bai = aligner + "/{sample}.bam.bai"
                output:
                    bam_out = "filtered_bam/{sample}.filtered.bam",
                    bai_out = "filtered_bam/{sample}.filtered.bam.bai",
                shell: """
                    ln -s ../{input.bam} {output.bam_out};
                    ln -s ../{input.bai} {output.bai_out}
                """


        rule sambamba_flagstat:
           input:
               aligner + "/{sample}.bam"
           output:
               "Sambamba/{sample}.markdup.txt"
           conda: CONDA_SAMBAMBA_ENV
           shell: """
               sambamba flagstat -p {input} > {output}
               """
