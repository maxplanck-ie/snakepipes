

### Filter the genome fasta for selected chromsomes for picard if the mode is allele-specific
rule filterFasta:
    input:
        bam = expand(mapping_prg+"/{sample}.bam", sample = samples[0]),
        fasta = genome_fasta
    output:
        temp(mapping_prg+"/genome_selectedChrs.fa")
    run:
        shell("samtools view -H {input.bam}"
              " | awk '$1 == \"@SQ\" {{print $2}}' | sed 's/SN://g' | sort -n > chrnames.txt")
        shell("samtools faidx {input.fasta} `cat chrnames.txt ` > {output} && rm chrnames.txt")


### Picard CollectAlignmentSummaryMetrics ######################################
<<<<<<< HEAD
## use filtered genome fasta if the mode is allele-specific

rule CollectAlignmentSummaryMetrics:
    input:
        bam = mapping_prg+"/{sample}.bam",
        genome =  lambda wildcards: mapping_prg+"/genome_selectedChrs.fa" if 'allelic-mapping' in mode else genome_fasta
    output:
        "Picard_qc/AlignmentSummaryMetrics/{sample}.alignment_summary_metrics.txt"
    #params:
        #genome =  lambda wildcards: mapping_prg+"/genome_selectedChrs.fa" if 'allelic-mapping' in mode else genome_fasta  # reference genome FASTA sequence
=======
## skip providing genome fasta if the mode is allele-specific

rule CollectAlignmentSummaryMetrics:
    input:
        mapping_prg+"/{sample}.bam"
    output:
        "Picard_qc/AlignmentSummaryMetrics/{sample}.alignment_summary_metrics.txt"
    params:
        genome =  lambda wildcards: 'null' if 'allelic-mapping' in mode else genome_fasta  # reference genome FASTA sequence
>>>>>>> 0d031db3eafd6a3ae5d5ba3c7d9242dbf3c28ef6
    log:
        "Picard_qc/logs/CollectAlignmentSummaryMetrics.{sample}.log"
    benchmark:
        "Picard_qc/.benchmark/CollectAlignmentSummaryMetrics.{sample}.benchmark"
    threads: 4 # Java performs parallel garbage collection
    shell:
        "java -Xmx4g -jar "+picard_path+"picard.jar CollectAlignmentSummaryMetrics "
            "REFERENCE_SEQUENCE={input.genome} "
            "INPUT={input.bam} OUTPUT={output} "
            "VALIDATION_STRINGENCY=LENIENT "
            "&> {log}"


### Picard CollectInsertSizeMetrics ############################################

if paired:
    rule CollectInsertSizeMetrics:
        input:
            mapping_prg+"/{sample}.bam"
        output:
            txt = "Picard_qc/InsertSizeMetrics/{sample}.insert_size_metrics.txt"
        log:
            "Picard_qc/logs/CollectInsertSizeMetrics.{sample}.log"
        benchmark:
            "Picard_qc/.benchmark/CollectInsertSizeMetrics.{sample}.benchmark"
        threads: 4 # Java performs parallel garbage collection
        params:
            pdf = "Picard_qc/InsertSizeMetrics/{sample}.insert_size_histogram.pdf"
        shell:
            "export PATH="+R_path+":$PATH && "
            "java -Xmx4g -jar "+picard_path+"picard.jar CollectInsertSizeMetrics "
                "HISTOGRAM_FILE={params.pdf} "
                "INPUT={input} "
                "OUTPUT={output.txt} "
                "VALIDATION_STRINGENCY=LENIENT "
                "&> {log} "
                "&& ( [ -f {params.pdf} ] || "+os.path.join(R_path, "Rscript")+" "+os.path.join(maindir, "shared", "tools", "CollectInsertSizeMetrics_histogram.R")+" {output.txt} ) "
