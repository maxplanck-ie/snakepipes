### samtools_filter ############################################################


# When modifying the rule samtools_filter, double-check whether the function
# update_filter() has to be modified too

rule samtools_filter:
    input:
        aligner+"/{sample}.bam",
        "filter_rules"
    output:
        bam = temp("filtered_bam/{sample}.filtered.tmp.bam")
    params:
        dedup = dedup,
        properPairs = properPairs,
        mapq = mapq
    log:
        out = "filtered_bam/logs/samtools_filter.{sample}.out",
        err = "filtered_bam/logs/samtools_filter.{sample}.err"
    benchmark:
        "filtered_bam/.benchmark/samtools_filter.{sample}.benchmark"
    threads: 8
    conda: CONDA_SHARED_ENV
    shell: """
        filter=""
        if [ "{params.dedup}" == "True" ] ; then filter="$filter -F 1024"; fi
        if [ "{params.properPairs}" == "True" ] ; then filter="$filter -f 2"; fi
        if [ "{params.mapq}" != "0" ] ; then filter="$filter -q {params.mapq}"; fi
        if [[ -z $filter ]] ; then ln -s -r {input[0]} {output.bam} ;
        else
            samtools view -@ {threads} -b $filter -o {output.bam} {input[0]} 2> {log.err} ;
        fi
        echo "samtools view arguments: $filter" > {log.out}
        """


### samtools_index #############################################################
rule samtools_index_tmp_filtered:
    input:
        "filtered_bam/{sample}.filtered.tmp.bam"
    output:
        temp("filtered_bam/{sample}.filtered.tmp.bam.bai")
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"
