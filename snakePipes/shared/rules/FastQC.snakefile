if pipeline == "scrnaseq" and mode == "STARsolo" or pipeline=="scrnaseq" and mode == "Alevin":
    rule FastQC:
        input:
            "originalFASTQ/{sample}{read}.fastq.gz"
        output:
            "FastQC/{sample}{read}_fastqc.html"
        benchmark:
            "FastQC/.benchmark/FastQC.{sample}{read}.benchmark"
        threads: 2
        conda: CONDA_SHARED_ENV
        shell: "fastqc -o FastQC {input}"

else:
    if pairedEnd:
        rule FastQC:
            input:
                "EXTERNAL_BAM/{sample}."+bamExt if fromBAM else "FASTQ/{sample}{read}.fastq.gz"
            output:
                "FastQC/{sample}_fastqc.html" if fromBAM else "FastQC/{sample}{read}_fastqc.html"
            benchmark:
                "FastQC/.benchmark/FastQC.{sample}{read}.benchmark"
            threads: 2
            conda: CONDA_SHARED_ENV
            shell: "fastqc -o FastQC {input}"

    else:
        rule FastQC_singleEnd:
            input:
                "EXTERNAL_BAM/{sample}."+bamExt if fromBAM else "FASTQ/{sample}"+reads[0]+".fastq.gz"
            output:
                "FastQC/{sample}_fastqc.html" if fromBAM else "FastQC/{sample}"+reads[0]+"_fastqc.html"
            params:
                reads=reads[0]
            benchmark:
                "FastQC/.benchmark/FastQC.{sample}"+reads[0]+".benchmark"
            threads: 2
            conda: CONDA_SHARED_ENV
            shell: """
                fastqc -o FastQC {input}
                """
