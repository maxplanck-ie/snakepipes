### samtools_filter ############################################################
CONDA_SHARED_ENV = "envs/shared_environment.yaml"

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
    threads: 8
    conda: CONDA_SHARED_ENV
    shell: """
        filter=""
        if [ "{params.dedup}" == "True" ] ; then filter="$filter -F 1024"; fi
        if [ "{params.properpairs}" == "True" ] ; then filter="$filter -f 2"; fi
        if [ "{params.mapq}" != "0" ] ; then filter="$filter -q params.mapq"; fi
        if [ "filter" == ""] ; then
            ln -s -r {input} {output.bam} ;
        else
            samtools view -@ {threads} -b $filter -q {{input} > {output.bam} 2> {log} ;
        fi
        echo 'samtools view arguments: $filter' > {output.filter_file}
        """


### samtools_index #############################################################
rule samtools_index_filtered:
    input:
        "filtered_bam/{sample}.bam"
    output:
        "filtered_bam/{sample}.bam.bai"
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"
