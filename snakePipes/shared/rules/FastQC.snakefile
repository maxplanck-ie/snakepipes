if pipeline == "scrnaseq" and mode == "STARsolo" or pipeline=="scrnaseq" and mode == "Alevin":
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
        conda: CONDA_SHARED_ENV
        shell: "fastqc -o FastQC {input} > {log.out} 2> {log.err}"

else:
    if pairedEnd:
        rule FastQC:
            input:
                "EXTERNAL_BAM/{sample}."+bamExt if fromBAM else "FASTQ/{sample}{read}.fastq.gz"
            output:
                "FastQC/{sample}_fastqc.html" if fromBAM else "FastQC/{sample}{read}_fastqc.html"
            log:
                out = "FastQC/logs/FastQC.{sample}{read}.out",
                err = "FastQC/logs/FastQC.{sample}{read}.err"
            benchmark:
                "FastQC/.benchmark/FastQC.{sample}{read}.benchmark"
            threads: 2
            conda: CONDA_SHARED_ENV
            shell: "fastqc -o FastQC {input} > {log.out} 2> {log.err}"

    else:
        rule FastQC_singleEnd:
            input:
                "EXTERNAL_BAM/{sample}."+bamExt if fromBAM else "FASTQ/{sample}"+reads[0]+".fastq.gz"
            output:
                "FastQC/{sample}_fastqc.html" if fromBAM else "FastQC/{sample}"+reads[0]+"_fastqc.html"
            params:
                reads=reads[0]
            log:
                out = "FastQC/logs/FastQC.{sample}"+reads[0]+".out",
                err = "FastQC/logs/FastQC.{sample}"+reads[0]+".err"
            benchmark:
                "FastQC/.benchmark/FastQC.{sample}"+reads[0]+".benchmark"
            threads: 2
            conda: CONDA_SHARED_ENV
            shell: """
                fastqc -o FastQC {input} > {log.out} 2> {log.err}
                """
