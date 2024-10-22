if pairedEnd or pipeline=="scrnaseq":
    rule validateFQ:
        input:
            r1 = indir+"/{sample}"+reads[0]+ext,
            r2 = indir+"/{sample}"+reads[1]+ext
        output:
            temp("originalFASTQ/{sample}.valid")
        conda: CONDA_FQLINT_ENV
        shell:"""
            fq lint {input.r1} {input.r2}
            touch {output}
            """
else:
    rule validateFQ:
        input:
            r1 = indir+"/{sample}"+reads[0]+ext
        output:
            temp("originalFASTQ/{sample}.valid")
        conda: CONDA_FQLINT_ENV
        shell:"""
            fq lint {input.r1}
            touch {output}
            """

rule origFASTQ1:
    input:
        r1 = indir+"/{sample}"+reads[0]+ext,
        valid = 'originalFASTQ/{sample}.valid'
    output:
        "originalFASTQ/{sample}"+reads[0]+".fastq.gz"
    params:
          cmd = lambda wildcards, input,output: "ln -s ../{} {}".format(input[0],output[0]) if pipeline=="preprocessing" else "ln -s {} {}".format(input[0],output[0])
    shell: """
               {params.cmd}
          """

if pairedEnd or pipeline=="scrnaseq":
    rule origFASTQ2:
        input:
            r2 = indir+"/{sample}"+reads[1]+ext,
            valid = 'originalFASTQ/{sample}.valid'
        output:
            "originalFASTQ/{sample}"+reads[1]+".fastq.gz"
        params:
            cmd = lambda wildcards, input,output: "ln -s ../{} {}".format(input[0],output[0]) if pipeline=="preprocessing" else "ln -s {} {}".format(input[0],output[0])
        shell: """
               {params.cmd}
            """


if downsample:
    if pairedEnd:
        rule FASTQdownsample:
            input:
                r1 = "originalFASTQ/{sample}"+reads[0]+".fastq.gz",
                r2 = "originalFASTQ/{sample}"+reads[1]+".fastq.gz"
            output:
                r1 = "originalFASTQ/downsample_{sample}"+reads[0]+".fastq.gz",
                r2 = "originalFASTQ/downsample_{sample}"+reads[1]+".fastq.gz"
            params:
                num_reads = downsample
            benchmark:
                "FASTQ/.benchmark/FASTQ_downsample.{sample}.benchmark"
            threads: lambda wildcards: 10 if 10<max_thread else max_thread
            conda: CONDA_SHARED_ENV
            shell: """
                seqtk sample -s 100 {input.r1} {params.num_reads} | pigz -p {threads} -9 > {output.r1}
                seqtk sample -s 100 {input.r2} {params.num_reads} | pigz -p {threads} -9 > {output.r2}
                """
    else:
        rule FASTQdownsample:
            input:
                "originalFASTQ/{sample}.fastq.gz"
            output:
                fq = "originalFASTQ/downsample_{sample}.fastq.gz"
            threads: lambda wildcards: 12 if 12<max_thread else max_thread
            params:
                num_reads = downsample
            conda: CONDA_SHARED_ENV
            shell: """
                seqtk sample -s 100 {input} {params.num_reads} | pigz -p {threads} -9 > {output}
                """
