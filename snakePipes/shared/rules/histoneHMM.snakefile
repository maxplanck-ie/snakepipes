### histoneHMM broad enrichment calling ########################################

def format_HMM_output(infile,outfile):
    f=open(infile)
    header=''
    body=[]
    for idx, line in enumerate(f):
        if idx==0:
            header=line
        if idx>0:
            cols = line.split("\t")
            cols1_3=cols[1:3]
            cols1_3_new=['{0:n}'.format(int(float(x))) for x in cols1_3]
            cols_new=cols
            cols_new[1:3]=cols1_3_new
            j='\t'.join(cols)
            body.append(j)
    f.close()
    with open(outfile, 'w') as f:
        f.write("%s" % header)
    with open(outfile, 'a') as f:
        for item in body:
            f.write("%s" % item)


# -b 750 -P 0.1
rule histoneHMM:
    input:
        "filtered_bam/{sample}.filtered.bam" if not useSpikeInForNorm else "split_bam/{sample}_host.bam"
    output:
        temp("histoneHMM/{sample}.filtered.histoneHMM-regions.gff"),
        temp("histoneHMM/{sample}.filtered.histoneHMM-em-posterior.txt"),
        temp("histoneHMM/{sample}.filtered.histoneHMM.txt")
    params:
        prefix = "histoneHMM/{sample}.filtered.histoneHMM",
        genome_index = genome_index
    benchmark:
        "histoneHMM/.benchmark/histoneHMM.{sample}.filtered.benchmark"
    conda: CONDA_CHIPSEQ_ENV
    shell: """
        RHOME=`R RHOME`
        $RHOME/library/histoneHMM/bin/histoneHMM_call_regions.R -b 750 -c {params.genome_index} -o {params.prefix} -P 0.1 {input}
        """

rule format_HMM_output:
    input:
        post ="histoneHMM/{sample}.filtered.histoneHMM-em-posterior.txt",
        txt = "histoneHMM/{sample}.filtered.histoneHMM.txt"
    output:
        post = temp("histoneHMM/{sample}.filtered.histoneHMM-em-posterior_formatted.txt"),
        txt = temp("histoneHMM/{sample}.filtered.histoneHMM_formatted.txt")
    run:
        format_HMM_output(input.post,output.post)
        format_HMM_output(input.txt,output.txt)


### compress and index GFF result file from histoneHMM for usage with IGV ######
### compress txt result files to save space ####################################
rule histoneHMM_out_gz:
    input:
        gff = "histoneHMM/{sample}.filtered.histoneHMM-regions.gff",
        post = "histoneHMM/{sample}.filtered.histoneHMM-em-posterior_formatted.txt",
        txt = "histoneHMM/{sample}.filtered.histoneHMM_formatted.txt"
    output:
        gff = "histoneHMM/{sample}.filtered.histoneHMM-regions.gff.gz",
        # touch output files as their modification date must more recent than
        # the modification date of the input files
        post = touch("histoneHMM/{sample}.filtered.histoneHMM-em-posterior.txt.gz"),
        txt = touch("histoneHMM/{sample}.filtered.histoneHMM.txt.gz")
    benchmark:
        "histoneHMM/.benchmark/histoneHMM_out_gz.{sample}.filtered.benchmark"
    threads: 2
    conda: CONDA_SHARED_ENV
    shell: """
        grep -v ^\"#\" {input.gff} | sort -k1,1 -k4,4n | bgzip > {output.gff}
        tabix -p gff {output.gff}
        gzip {input.post}
        gzip {input.txt}
        """
