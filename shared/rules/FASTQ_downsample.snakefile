rule fastq_downlsample:
    input:  indir+"/{sample}{read}.fastq.gz"
    output: "FASTQ/{sample}{read}.fastq.gz"
    benchmark:  "FASTQ/.benchmark/FASTQ_downsample.{sample}{read}.benchmark"
    threads: 4
    shell:  "bash %s {threads} %i {input} {output}" % (os.path.join(workflow_tools,"fastq_head.sh"), int(config["downsample"]))
