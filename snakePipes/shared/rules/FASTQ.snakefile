import os
rule origFASTQ1:
      input:
          indir+"/{sample}"+reads[0]+ext
      output:
          "originalFASTQ/{sample}"+reads[0]+".fastq.gz"
      run:
        if not os.path.exists(os.path.join(outdir,output[0])):
            os.symlink(os.path.join(outdir,input[0]),os.path.join(outdir,output[0]))

if pairedEnd or pipeline=="scRNA-seq":
    rule origFASTQ2:
        input:
            indir+"/{sample}"+reads[1]+ext
        output:
            "originalFASTQ/{sample}"+reads[1]+".fastq.gz"
        run:
            if not os.path.exists(os.path.join(outdir,output[0])):
                os.symlink(os.path.join(outdir,input[0]),os.path.join(outdir,output[0]))

if downsample:
    if pairedEnd:
        rule FASTQdownsample:
            input:
                r1 = "originalFASTQ/{sample}"+reads[0]+".fastq.gz",
                r2 = "originalFASTQ/{sample}"+reads[1]+".fastq.gz"
            output:
                r1 = "originalFASTQ/downsample_{sample}"+reads[0]+".fastq.gz",
                r2 = "originalFASTQ/downsample_{sample}"+reads[1]+".fastq.gz"
            log: "originalFASTQ/logs/{sample}.FASTQdownsample.log"
            params:
                num_reads = downsample
            benchmark:
                "FASTQ/.benchmark/FASTQ_downsample.{sample}.benchmark"
            threads: 10
            conda: CONDA_SHARED_ENV
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
            threads: 12
            params:
                num_reads = downsample
            conda: CONDA_SHARED_ENV
            shell: """
                seqtk sample -s 100 {input} {params.num_reads} | pigz -p {threads} -9 > {output} 2> {log}
                """
