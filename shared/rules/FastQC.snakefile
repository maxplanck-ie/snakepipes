rule FastQC:
    input:  "FASTQ/{sample}{read}.fastq.gz"
    output: "FastQC/{sample}{read}_fastqc.html"
    log:    "FastQC/log/{sample}{read}.log"
    benchmark:  "FastQC/.benchmark/FastQC.{sample}{read}.benchmark"
    shell:  fastqc_path+" -o FastQC {input} 2>&1 | tee {log}"
