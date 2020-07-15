import subprocess

part=['host','spikein']

# MACS2 should be called on already filtered, e.g. duplicate-free, BAM files
# for paired-end BAM files, sambamba markdupes is fragment-based and
# therefore superior to MACS2 mate 1-based duplicate detection


### MACS2 peak calling #########################################################

if pairedEnd:
    rule writeFragmentSize:
        input: "split_deepTools_qc/bamPEFragmentSize/host.fragmentSize.metric.tsv"
        output: "MACS2/fragmentSize.metrix.tsv"


    rule MACS2:
        input:
            chip = "split_bam/{chip_sample}_host.bam",
            control =
                lambda wildcards: "split_bam/"+get_control(wildcards.chip_sample)+"_host.bam" if get_control(wildcards.chip_sample)
                else [],
            insert_size_metrics = "split_deepTools_qc/bamPEFragmentSize/host.fragmentSize.metric.tsv"
        output:
            peaks = "MACS2/{chip_sample}_host.BAM_peaks.xls",
            peaksPE = "MACS2/{chip_sample}_host.BAMPE_peaks.xls"
        params:
            genome_size = genome_size,
            broad_calling =
                lambda wildcards: "--broad" if is_broad(wildcards.chip_sample) else "",
            control_param =
                lambda wildcards: "-c split_bam/"+get_control(wildcards.chip_sample)+"_host.bam" if get_control(wildcards.chip_sample)
                else "",
            qval_cutoff=qval,
            mfold=mfold
        log:
            out = "MACS2/logs/MACS2.{chip_sample}_host.filtered.out",
            err = "MACS2/logs/MACS2.{chip_sample}_host.filtered.err"
        benchmark:
            "MACS2/.benchmark/MACS2.{chip_sample}_host.filtered.benchmark"
        conda: CONDA_CHIPSEQ_ENV
        shell: """
            macs2 callpeak -t {input.chip} {params.control_param} \
                -f BAM \
                -g {params.genome_size} --qvalue {params.qval_cutoff}\
                --keep-dup all \
                --outdir MACS2 \
                --name {wildcards.chip_sample}_host.BAM \
                --nomodel \
                --mfold {params.mfold}\
                --extsize $(cat {input.insert_size_metrics} | grep split_bam/{wildcards.chip_sample}_host.bam | awk '{{printf("%i",$6)}}') \
                {params.broad_calling} > {log.out} 2> {log.err}

            # also run MACS2 in paired-end mode BAMPE for comparison with single-end mode
            macs2 callpeak -t {input.chip} \
                {params.control_param} -f BAMPE --qvalue {params.qval_cutoff}\
                -g {params.genome_size} --keep-dup all \
                --outdir MACS2 --name {wildcards.chip_sample}_host.BAMPE \
                {params.broad_calling} > {log.out}.BAMPE 2> {log.err}.BAMPE
            """
else:
    rule MACS2:
        input:
            chip = "split_bam/{chip_sample}_host.bam",
            control =
                lambda wildcards: "split_bam/"+get_control(wildcards.chip_sample)+"_host.bam" if get_control(wildcards.chip_sample)
                else []
        output:
            peaks = "MACS2/{chip_sample}_host.BAM_peaks.xls",
        params:
            genome_size = int(genome_size),
            broad_calling =
                lambda wildcards: "--broad" if is_broad(wildcards.chip_sample)
                else "",
            control_param =
                lambda wildcards: "-c split_bam/"+get_control(wildcards.chip_sample)+"_host.bam" if get_control(wildcards.chip_sample)
                else "",
            frag_size=fragmentLength,
            mfold=mfold,
            qval_cutoff=qval
        log:
            out = "MACS2/logs/MACS2.{chip_sample}_host.filtered.out",
            err = "MACS2/logs/MACS2.{chip_sample}_host.filtered.err"
        benchmark:
            "MACS2/.benchmark/MACS2.{chip_sample}_host.filtered.benchmark"
        conda: CONDA_CHIPSEQ_ENV
        shell: """
            macs2 callpeak -t {input.chip} {params.control_param} -f BAM -g {params.genome_size} --qvalue {params.qval_cutoff} --keep-dup all --outdir MACS2 \
                --name {wildcards.chip_sample}_host.BAM --mfold {params.mfold} --extsize {params.frag_size}\
                {params.broad_calling} > {log.out} 2> {log.err}
            """



### MACS2 peak quality control #################################################

rule MACS2_peak_qc:
    input:
        bam = "split_bam/{sample}_host.bam",
        xls = "MACS2/{sample}_host.bam_peaks.xls"
    output:
        qc = "MACS2/{sample}_host.bam_peaks.qc.txt"
    params:
        peaks =
            lambda wildcards: "MACS2/{}_host.BAM_peaks.broadPeak".format(wildcards.sample) if is_broad(wildcards.sample)
                              else "MACS2/{}_host.BAM_peaks.narrowPeak".format(wildcards.sample),
        genome_index = genome_index
    benchmark:
        "MACS2/.benchmark/MACS2_peak_qc.{sample}_host.filtered.benchmark"
    conda: CONDA_SHARED_ENV
    shell: """
        # get the number of peaks
        peak_count=`wc -l < {params.peaks}`

        # get the number of mapped reads
        mapped_reads=`samtools view -c -F 4 {input.bam}`

        # calculate the number of alignments overlapping the peaks
        # exclude reads flagged as unmapped (unmapped reads will be reported when using -L)
        reads_in_peaks=`samtools view -c -F 4 -L {params.peaks} {input.bam}`

        # calculate Fraction of Reads In Peaks
        frip=`bc -l <<< "$reads_in_peaks/$mapped_reads"`

        # compute peak genome coverage
        peak_len=`awk '{{total+=$3-$2}}END{{print total}}' {params.peaks}`
        genome_size=`awk '{{total+=$3-$2}}END{{print total}}' {params.genome_index}`
        genomecov=`bc -l <<< "$peak_len/$genome_size"`

        # write peak-based QC metrics to output file
        printf "peak_count\tFRiP\tpeak_genome_coverage\n%d\t%5.3f\t%6.4f\n" $peak_count $frip $genomecov > {output.qc}
        """

# TODO
# add joined deepTools plotEnrichment call for all peaks and samples in one plot

rule namesort_bams:
    input:
        bam = "split_bam/{sample}.host.bam"
    output:
        bam = temp("split_bam/{sample}.namesorted.bam")
    log:
        "split_bam/logs/namesort.err"
    params:
        tempDir = tempDir
    threads: 4
    conda: CONDA_SAMBAMBA_ENV
    shell: """
        TMPDIR={params.tempDir}
        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX)
        sambamba sort -t {threads} -o {output.bam} --tmpdir=$MYTEMP -n {input.bam} 2> {log}
        rm -rf $MYTEMP
         """

# Requires PE data
# Should be run once per-group!
if pairedEnd:
    rule Genrich_peaks:
        input:
            bams=lambda wildcards: expand(os.path.join("split_bam", "{sample}_namesorted.bam"), sample=genrichDict[wildcards.group]),
            control = lambda wildcards: ["split_bam/"+get_control(x)+"_namesorted.bam" for x in genrichDict[wildcards.group]]
        output:
            "Genrich/{group}.narrowPeak"
        log: "Genrich/logs/{group}.log"
        params:
            bams = lambda wildcards: ",".join(expand(os.path.join("split_bam", "{sample}_namesorted.bam"), sample=genrichDict[wildcards.group])),
            blacklist = "-E {}".format(blacklist_bed) if blacklist_bed else "",
            control_pfx=lambda wildcards,input: "-c" if input.control else "",
            control=lambda wildcards,input: ",".join(input.control) if input.control else "",
            spikein_chroms=",".join(spikein_chr)
        conda: CONDA_CHIPSEQ_ENV
        shell: """
            Genrich -S -t {params.bams} {params.control_pfx} {params.control} -o {output} -r {params.blacklist} -e {params.spikein_chroms} -y 2> {log}
            """
else:
    rule Genrich_peaks:
        input:
            bams=lambda wildcards: expand(os.path.join("split_bam", "{sample}_namesorted.bam"), sample=genrichDict[wildcards.group]),
            control = lambda wildcards: ["split_bam/"+get_control(x)+"_namesorted.bam" for x in genrichDict[wildcards.group]]
        output:
            "Genrich/{group}.narrowPeak"
        log: "Genrich/logs/{group}.log"
        params:
            bams = lambda wildcards: ",".join(expand(os.path.join("split_bam", "{sample}_namesorted.bam"), sample=genrichDict[wildcards.group])),
            blacklist = "-E {}".format(blacklist_bed) if blacklist_bed else "",
            control_pfx=lambda wildcards,input: "-c" if input.control else "",
            control=lambda wildcards,input: ",".join(input.control) if input.control else "",
            frag_size=fragmentLength,
            spikein_chroms=",".join(spikein_chr)
        conda: CONDA_CHIPSEQ_ENV
        shell: """
            Genrich -S -t {params.bams} {params.control_pfx} {params.control} -o {output} -r {params.blacklist} -e {params.spikein_chroms} -w {params.frag_size} 2> {log}
            """
