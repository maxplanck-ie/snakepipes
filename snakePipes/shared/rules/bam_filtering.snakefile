### samtools_filter ############################################################

# When modifying the rule samtools_filter, double-check whether the function
# update_filter() has to be modified too

bam_filter_string = "{} {} {}".format("-F 1024" if dedup else "", "-f 2" if properPairs else "", "-q "+ str(mapq) if str(mapq) != "0" else "")

rule samtools_filter:
    input:
        aligner+"/{sample}.bam"
    output:
        bam = temp("filtered_bam/{sample}.filtered.tmp.bam")
    params:
        shell = lambda wildcards,input,output: "samtools view -@ {} -b {} -o {} {} ".format(str(8 if not local else 2), bam_filter_string,output.bam,input[0]) if bam_filter_string.strip() !="" else "ln -s ../{} {}".format(input[0],output.bam)
    benchmark:
        "filtered_bam/.benchmark/samtools_filter.{sample}.benchmark"
    threads: lambda wildcards: 8 if 8<max_thread else max_thread
    conda: CONDA_SHARED_ENV
    shell: """
        {params.shell}
        """


### samtools_index #############################################################
rule samtools_index_tmp_filtered:
    input:
        "filtered_bam/{sample}.filtered.tmp.bam"
    output:
        temp("filtered_bam/{sample}.filtered.tmp.bam.bai")
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"
