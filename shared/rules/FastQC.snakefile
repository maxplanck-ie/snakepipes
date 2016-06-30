rule FastQC:
    input:
        "FASTQ/{sample}{read}.fastq.gz"
    output:
        "FastQC/{sample}{read}_fastqc.html"
    log:
        "FastQC/logs/FastQC.{sample}{read}.log"
    benchmark:
        "FastQC/.benchmark/FastQC.{sample}{read}.benchmark"
    threads: 2
    shell:
        fastqc_path+"fastqc -o FastQC {input} &> {log}"
