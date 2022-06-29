import os

rule link_bam:
    input:
        config['indir'] + "/{sample}" + config['bamExt']
    output:
        config['aligner'] + "/{sample}.unsorted.bam" if config['pipeline']=="noncoding-rna-seq" else config['aligner'] + "/{sample}.bam"
    params:
        input_bai = config['indir'] + "/{sample}" + config['bamExt'] + ".bai",
        output_bai = config['aligner'] + "/{sample}.unsorted.bam.bai" if config['pipeline']=="noncoding-rna-seq" else config['aligner'] + "/{sample}.bam.bai"
    run:
        if os.path.exists(params.input_bai) and not os.path.exists(os.path.join(outdir,params.output_bai)):
            os.symlink(params.input_bai,os.path.join(outdir,params.output_bai))
        if not os.path.exists(os.path.join(outdir,output[0])):
            os.symlink(os.path.join(outdir,input[0]),os.path.join(outdir,output[0]))

if not config['pipeline']=="noncoding-rna-seq":
    rule samtools_index_external:
        input:
            config['aligner'] + "/{sample}.bam"
        output:
            config['aligner'] + "/{sample}.bam.bai"
        conda: config['CONDA_SHARED_ENV']
        shell: "if [[ ! -f {output[0]} ]]; then samtools index {input[0]}; fi"

    if not config['pipeline']=="WGBS" or config['pipeline']=="WGBS" and skipBamQC:
        rule link_bam_bai_external:
            input:
                bam = config['aligner'] + "/{sample}.bam",
                bai = config['aligner'] + "/{sample}.bam.bai"
            output:
                bam_out = "filtered_bam/{sample}.filtered.bam",
                bai_out = "filtered_bam/{sample}.filtered.bam.bai",
            shell: """
                ln -s ../{input.bam} {output.bam_out};
                ln -s ../{input.bai} {output.bai_out}
            """


    rule sambamba_flagstat:
       input:
           config['aligner'] + "/{sample}.bam"
       output:
           "Sambamba/{sample}.markdup.txt"
       log: "Sambamba/logs/{sample}.flagstat.log"
       conda: config['CONDA_SAMBAMBA_ENV']
       shell: """
           sambamba flagstat -p {input} > {output} 2> {log}
           """
