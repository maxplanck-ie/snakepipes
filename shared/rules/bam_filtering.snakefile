### Picard MarkDuplicates ######################################################

rule MarkDuplicates:
    input:
        mapping_prg+"/{sample}.sorted.bam"
    output:
        bam = mapping_prg+"/{sample}.bam",
        txt = "Picard_qc/MarkDuplicates/{sample}.mark_duplicates_metrics.txt"
    log:
        "Picard_qc/logs/MarkDuplicates.{sample}.log"
    benchmark:
        "Picard_qc/.benchmark/MarkDuplicates.{sample}.benchmark"
    threads: 4 # Java performs parallel garbage collection
    shell:
        "java -Xmx8g -jar "+picard_path+"picard.jar MarkDuplicates "
            "MAX_FILE_HANDLES=1000 "
            "INPUT={input} "
            "OUTPUT={output.bam} "
            "METRICS_FILE={output.txt} "
            "ASSUME_SORTED=true "
            "VERBOSITY=WARNING "
            "VALIDATION_STRINGENCY=LENIENT "
            "&> {log} "


### samtools_filter ############################################################
# When modifying the rule samtools_filter, double-check wether the function
# update_filter() has to be modified concordantly

rule samtools_filter:
    input:
        mapping_prg+"/{sample}.bam"
    output:
        bam = "filtered_bam/{sample}.filtered.bam",
        filter_file = "filtered_bam/{sample}.filter" if (dedup or properpairs or mapq > 0)
                      else []
    params:
        dedup = dedup,
        properpairs = properpairs,
        mapq = mapq
    log:
        "filtered_bam/logs/samtools_filter.{sample}.log"
    benchmark:
        "filtered_bam/.benchmark/samtools_filter.{sample}.benchmark"
    run:
        # string with samtools view parameters for filtering
        filter = ""
        if params.dedup:
            filter += "-F 1024 "
        if params.properpairs:
            filter += "-f 2 "
        if params.mapq > 0:
            filter += "-q {params.mapq}"

        if filter:
            shell(
                samtools_path+"samtools view "
                "-b "+filter+" {input} > {output.bam} "
                "2> {log} "
                "&& echo 'samtools view arguments: "+filter+"' > {output.filter_file}"
            )
        else:
            shell(
                "( [ -f {output.bam} ] || ln -s -r {input} {output.bam} ) && touch -h {output.bam}"
            )


### samtools_index #############################################################

rule samtools_index:
    input:
        mapping_prg+"/{sample}.bam"
    output:
        mapping_prg+"/{sample}.bam.bai"
    shell:
        samtools_path+"samtools index {input}"

rule samtools_index2:
    input:
        "filtered_bam/{sample}.bam"
    output:
        "filtered_bam/{sample}.bam.bai"
    shell:
        samtools_path+"samtools index {input}"
