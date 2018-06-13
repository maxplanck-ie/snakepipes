

if paired:
    rule FastQC:
        input:
            "FASTQ/{sample}{read}.fastq.gz"
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
    rule FastQC_singleEnd:
        input:
            "FASTQ/{sample}.fastq.gz"
        output:
            "FastQC/{sample}_fastqc.html"
        log:
            out = "FastQC/logs/FastQC.{sample}.out",
            err = "FastQC/logs/FastQC.{sample}.err"
        benchmark:
            "FastQC/.benchmark/FastQC.{sample}.benchmark"
        threads: 2
        conda: CONDA_SHARED_ENV
        shell: "fastqc -o FastQC {input} > {log.out} 2> {log.err}"
