## picard: CollectAlignmentSummaryMetrics ######################################

rule CollectAlignmentSummaryMetrics:
    input:
        fasta = genome_fasta,
        bam = "Bowtie2/{sample}.bam"
    output: "AlignmentSummaryMetrics/{sample}.AlignmentSummaryMetrics.txt"
    log:    "AlignmentSummaryMetrics/log/{sample}.CollectAlignmentSummaryMetrics.log"
    benchmark:  "AlignmentSummaryMetrics/.benchmark/Picard_CollectAlignmentSummaryMetrics.{sample}.benchmark"
    shell:  "java -Xmx4g -jar "+picard_path+" CollectAlignmentSummaryMetrics "
            "VALIDATION_STRINGENCY=LENIENT "
            "R={input.fasta} "
            "I={input.bam} O={output} "
            "2>&1 | tee {log} "


## picard: CollectInsertSizeMetrics ######################################

if paired:
    rule CollectInsertSizeMetrics:
        input:   "Bowtie2/{sample}.bam"
        output:
            txt = "InsertSizeMetrics/{sample}.InsertSizeMetrics.txt",
            pdf = "InsertSizeMetrics/{sample}.InsertSizeMetrics.pdf",
        log:    "InsertSizeMetrics/log/{sample}.CollectInsertSizeMetrics.log"
        benchmark:  "InsertSizeMetrics/.benchmark/Picard_CollectInsertSizeMetrics.{sample}.benchmark"
        shell:  "export PATH="+R_dir+":$PATH && "
                "java -Xmx4g -jar "+picard_path+" CollectInsertSizeMetrics "
                "VALIDATION_STRINGENCY=LENIENT "
                "I={input} "
                "HISTOGRAM_FILE={output.pdf} O={output.txt} "
                "2>&1 | tee {log} "

    rule CollectInsertSizeMetrics_extract_mean:
        input:  "InsertSizeMetrics/{sample}.InsertSizeMetrics.txt"
        output: "InsertSizeMetrics/{sample}.mean.txt"
        shell:  """mean=$(cat {input} | head -8 | tail -1 | awk '{{printf "mean\t%.0f", $5}}'); echo $mean > {output}"""
