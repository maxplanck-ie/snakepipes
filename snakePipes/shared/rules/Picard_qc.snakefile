


### Filter the genome fasta for selected chromsomes for picard if the mode is allele-specific
rule filterFasta:
    input:
        bam = expand(mapping_prg+"/{sample}.bam", sample = samples[0]),
        fasta = genome_fasta
    output:
        temp(mapping_prg+"/genome_selectedChrs.fa")
    conda: CONDA_SHARED_ENV
    shell: """
        samtools view -H {input.bam} | awk '$1 == \"@SQ\" {{print $2}}' | sed 's/SN://g' | sort -n > chrnames.txt
        samtools faidx {input.fasta} `cat chrnames.txt ` > {output} && rm chrnames.txt
        """


### Picard CollectAlignmentSummaryMetrics ######################################

# use filtered genome fasta if the mode is allele-specific


rule CollectAlignmentSummaryMetrics:
    input:
        bam = mapping_prg+"/{sample}.bam",
        genome =  lambda wildcards: mapping_prg+"/genome_selectedChrs.fa" if 'allelic-mapping' in mode else genome_fasta
    output:
        "Picard_qc/AlignmentSummaryMetrics/{sample}.alignment_summary_metrics.txt"
    log:
        out = "Picard_qc/logs/CollectAlignmentSummaryMetrics.{sample}.out",
        err = "Picard_qc/logs/CollectAlignmentSummaryMetrics.{sample}.err"
    benchmark:
        "Picard_qc/.benchmark/CollectAlignmentSummaryMetrics.{sample}.benchmark"
    threads: 4 # Java performs parallel garbage collection, 2GB per core
    conda: CONDA_DNA_MAPPING_ENV
    shell:
        "picard CollectAlignmentSummaryMetrics "
        "REFERENCE_SEQUENCE={input.genome} "
        "INPUT={input.bam} OUTPUT={output} "
        "VALIDATION_STRINGENCY=LENIENT "
        "> {log.out} 2> {log.err}"


### Picard CollectInsertSizeMetrics ############################################

if paired:
    rule CollectInsertSizeMetrics:
        input:
            mapping_prg+"/{sample}.bam"
        output:
            txt = "Picard_qc/InsertSizeMetrics/{sample}.insert_size_metrics.txt"
        log:
            out = "Picard_qc/logs/CollectInsertSizeMetrics.{sample}.out",
            err = "Picard_qc/logs/CollectInsertSizeMetrics.{sample}.err"
        benchmark:
            "Picard_qc/.benchmark/CollectInsertSizeMetrics.{sample}.benchmark"
        threads: 4 # Java performs parallel garbage collection, 1GB per core
        params:
            pdf = "Picard_qc/InsertSizeMetrics/{sample}.insert_size_histogram.pdf"
        conda: CONDA_DNA_MAPPING_ENV
        shell:
            "picard CollectInsertSizeMetrics "
            "HISTOGRAM_FILE={params.pdf} "
            "INPUT={input} "
            "OUTPUT={output.txt} "
            "VALIDATION_STRINGENCY=LENIENT "
            "> {log.out} 2> {log.err} "
