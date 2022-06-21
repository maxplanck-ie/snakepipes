if config['pipeline'] == "scrna-seq" and config['mode'] == "STARsolo" or config['pipeline']=="scrna-seq" and config['mode'] == "Alevin":
    rule FastQC:
        input:
            "originalFASTQ/{sample}{read}.fastq.gz"
        output:
            "FastQC/{sample}{read}_fastqc.html"
        log:
            out = "FastQC/logs/FastQC.{sample}{read}.out",
            err = "FastQC/logs/FastQC.{sample}{read}.err"
        benchmark:
            "FastQC/.benchmark/FastQC.{sample}{read}.benchmark"
        threads: 2
        conda: config['CONDA_SHARED_ENV']
        shell: "fastqc -o FastQC {input} > {log.out} 2> {log.err}"

else:
    if config['pairedEnd']:
        rule FastQC:
            input:
                "EXTERNAL_BAM/{sample}."+config['bamExt'] if config['fromBAM'] else "FASTQ/{sample}{read}.fastq.gz"
            output:
                "FastQC/{sample}_fastqc.html" if config['fromBAM'] else "FastQC/{sample}{read}_fastqc.html"
            log:
                out = "FastQC/logs/FastQC.{sample}{read}.out",
                err = "FastQC/logs/FastQC.{sample}{read}.err"
            benchmark:
                "FastQC/.benchmark/FastQC.{sample}{read}.benchmark"
            threads: 2
            conda: config['CONDA_SHARED_ENV']
            shell: "fastqc -o FastQC {input} > {log.out} 2> {log.err}"

    else:
        rule FastQC_singleEnd:
            input:
                "EXTERNAL_BAM/{sample}."+config['bamExt'] if config['fromBAM'] else "FASTQ/{sample}"+config['reads'][0]+".fastq.gz"
            output:
                "FastQC/{sample}_fastqc.html" if config['fromBAM'] else "FastQC/{sample}"+config['reads'][0]+"_fastqc.html"
            params:
                reads=reads[0]
            log:
                out = "FastQC/logs/FastQC.{sample}"+config['reads'][0]+".out",
                err = "FastQC/logs/FastQC.{sample}"+config['reads'][0]+".err"
            benchmark:
                "FastQC/.benchmark/FastQC.{sample}"+config['reads'][0]+".benchmark"
            threads: 2
            conda: config['CONDA_SHARED_ENV']
            shell: """
                fastqc -o FastQC {input} > {log.out} 2> {log.err}
                """
