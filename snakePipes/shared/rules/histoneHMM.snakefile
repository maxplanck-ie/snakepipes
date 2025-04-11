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


rule cleanup_histoneHMM:
    input:
        peaks = "histoneHMM/{sample}.filtered.histoneHMM-regions.gff"
    output:
        peaks_gff = "histoneHMM/{sample}_avgp0.5.gff",
        peaks_bed = "histoneHMM/{sample}_avgp0.5.bed"
    params:
        outdir = "histoneHMM",
        input_peaks = "../histoneHMM/{sample}.filtered.histoneHMM-regions.gff"
    conda: CONDA_CHIPQC_ENV
    script: "../rscripts/clean_histoneHMM_result.R"


rule histoneHMM_chipqc:
    input:
        bams = expand("filtered_bam/{broad_sample}.filtered.bam",broad_sample=broad_samples),
        peaks = expand("histoneHMM/{broad_sample}_avgp0.5.bed",broad_sample=broad_samples),
        sampleSheet = sampleSheet if sampleSheet else [],
        chipdict = os.path.join(outdir,"chip_samples.yaml")
    output:
        "histoneHMM_chipqc/sessionInfo.txt"
    params:
        genome = genome,
        outdir = "histoneHMM_chipqc",
        blacklist = blacklist_bed,
        bams = lambda wildcards,input: [os.path.join(outdir,x) for x in input.bams],
        peaks = lambda wildcards,input: [os.path.join(outdir,x) for x in input.peaks],
        narrow_samples = [],
        broad_samples = broad_samples,
        useSpikeinForNorm = useSpikeInForNorm
    threads: 8
    benchmark:
        "histoneHMM_chipqc/.benchmark/chipqc.benchmark"
    conda: CONDA_CHIPQC_ENV
    script: "../rscripts/chipqc.R"


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
        grep -v ^\"#\" {input.gff} | sort -k1,1 -k4,4n | bgzip -f  > {output.gff}
        tabix -f -p gff {output.gff}
        gzip -f {input.post}
        gzip -f {input.txt}
        """
