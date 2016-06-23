### histoneHMM broad enrichment calling ########################################

rule histoneHMM:
    input:
        "filtered_bam/{sample}.filtered.bam"
    output:
        temp("histoneHMM/{sample}.filtered.histoneHMM-regions.gff"),
        "histoneHMM/{sample}.filtered.histoneHMM-em-posterior.txt",
        "histoneHMM/{sample}.filtered.histoneHMM.txt"
    params:
        prefix = "histoneHMM/{sample}.filtered.histoneHMM"
    log:
        "histoneHMM/logs/histoneHMM.{sample}.filtered.log"
    benchmark:
        "histoneHMM/.benchmark/histoneHMM.{sample}.filtered.benchmark"
    shell:
      "export R_LIBS_USER="+R_libs_path+" &&"+
      R_path+"Rscript "+histoneHMM_path+"histoneHMM_call_regions.R "
        "-b 750 "
        "-c "+genome_index+" "
        "-o {params.prefix} "
        "-P 0.1 "
        "{input} "
        "&> {log}"


### compress and index GFF result file from histoneHMM for usage with IGV
### compress txt result files to save space
rule histoneHMM_out_gz:
    input:
        gff = "histoneHMM/{sample}.filtered.histoneHMM-regions.gff",
        post = "histoneHMM/{sample}.filtered.histoneHMM-em-posterior.txt",
        txt = "histoneHMM/{sample}.filtered.histoneHMM.txt"
    output:
        gff = "histoneHMM/{sample}.filtered.histoneHMM-regions.gff.gz",
        # touch output files as their modification date must more recent than
        # the modification date of the input files
        post = touch("histoneHMM/{sample}.filtered.histoneHMM-em-posterior.txt.gz"),
        txt = touch("histoneHMM/{sample}.filtered.histoneHMM.txt.gz")
    log:
        "histoneHMM/logs/histoneHMM_out_gz.{sample}.filtered.log"
    benchmark:
        "histoneHMM/.benchmark/histoneHMM_out_gz.{sample}.filtered.benchmark"
    threads: 2
    shell:
        "grep -v ^\"#\" {input.gff} | "
        "sort -k1,1 -k4,4n | "+
        tabix_path+"bgzip > {output.gff} && "+
        tabix_path+"tabix -p gff {output.gff} "
        "&> {log} && "+
        "gzip {input.post} && "+
        "gzip {input.txt}"
