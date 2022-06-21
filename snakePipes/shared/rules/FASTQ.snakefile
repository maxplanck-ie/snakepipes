rule origFASTQ1:
      input:
          config['indir']+"/{sample}"+config['reads'][0]+config['ext']
      output:
          "originalFASTQ/{sample}"+config['reads'][0]+".fastq.gz"
      params:
            cmd = lambda wildcards, input,output: "ln -s ../{} {}".format(input[0],output[0]) if config['pipeline']=="preprocessing" else "ln -s {} {}".format(input[0],output[0])
      shell: """
               {params.cmd}
          """

if config['pairedEnd'] or config['pipeline']=="scrna-seq":
    rule origFASTQ2:
        input:
            config['indir']+"/{sample}"+config['reads'][1]+config['ext']
        output:
            "originalFASTQ/{sample}"+config['reads'][1]+".fastq.gz"
        params:
            cmd = lambda wildcards, input,output: "ln -s ../{} {}".format(input[0],output[0]) if config['pipeline']=="preprocessing" else "ln -s {} {}".format(input[0],output[0])
        shell: """
               {params.cmd}
            """

if config['downsample']:
    if config['pairedEnd']:
        rule FASTQdownsample:
            input:
                r1 = "originalFASTQ/{sample}"+config['reads'][0]+".fastq.gz",
                r2 = "originalFASTQ/{sample}"+config['reads'][1]+".fastq.gz"
            output:
                r1 = "originalFASTQ/downsample_{sample}"+config['reads'][0]+".fastq.gz",
                r2 = "originalFASTQ/downsample_{sample}"+config['reads'][1]+".fastq.gz"
            log: "originalFASTQ/logs/{sample}.FASTQdownsample.log"
            params:
                num_reads = config['downsample']
            benchmark:
                "FASTQ/.benchmark/FASTQ_downsample.{sample}.benchmark"
            threads: lambda wildcards: 10 if 10<config['max_thread'] else config['max_thread']
            conda: config['CONDA_SHARED_ENV']
            shell: """
                seqtk sample -s 100 {input.r1} {params.num_reads} | pigz -p {threads} -9 > {output.r1} 2> {log}
                seqtk sample -s 100 {input.r2} {params.num_reads} | pigz -p {threads} -9 > {output.r2} 2>> {log}

                """
    else:
        rule FASTQdownsample:
            input:
                "originalFASTQ/{sample}.fastq.gz"
            output:
                fq = "originalFASTQ/downsample_{sample}.fastq.gz"
            log: "originalFASTQ/logs/{sample}.FASTQdownsample.log"
            threads: lambda wildcards: 12 if 12<config['max_thread'] else config['max_thread']
            params:
                num_reads = config['downsample']
            conda: config['CONDA_SHARED_ENV']
            shell: """
                seqtk sample -s 100 {input} {params.num_reads} | pigz -p {threads} -9 > {output} 2> {log}
                """
