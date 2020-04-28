### samtools_filter ############################################################
import os

# When modifying the rule samtools_filter, double-check whether the function
# update_filter() has to be modified too

bam_filter_string = "{} {} {}".format("-F 1024" if dedup else "", "-f 2" if properPairs else "", "-q "+ str(mapq) if str(mapq) != "0" else "")

rule samtools_filter:
    input:
        aligner+"/{sample}.bam"
    output:
        bam = temp("filtered_bam/{sample}.filtered.tmp.bam")
    params:
        shell = lambda wildcards,input,output: "samtools view -@ {} -b {} -o {} {} ".format(str(8 if not local else 2), bam_filter_string,output.bam,input[0]) if bam_filter_string.strip() !="" else "ln -s {} {}".format(os.path.join(outdir,input[0]),os.path.join(outdir,"filtered_bam",wildcards.sample+".filtered.tmp.bam")) 
    log:
        out = "filtered_bam/logs/samtools_filter.{sample}.out",
        err = "filtered_bam/logs/samtools_filter.{sample}.err"
    benchmark:
        "filtered_bam/.benchmark/samtools_filter.{sample}.benchmark"
    threads: 8
    conda: CONDA_SHARED_ENV
    shell: """
        {params.shell} 2> {log.err}
        """


### samtools_index #############################################################
rule samtools_index_tmp_filtered:
    input:
        "filtered_bam/{sample}.filtered.tmp.bam"
    output:
        temp("filtered_bam/{sample}.filtered.tmp.bam.bai")
    log: "filtered_bam/logs/{sample}.samtools_index_tmp_filtered.log"
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input} 2> {log}"
