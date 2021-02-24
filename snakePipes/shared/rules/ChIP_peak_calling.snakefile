import subprocess

# MACS2 should be called on already filtered, e.g. duplicate-free, BAM files
# for paired-end BAM files, sambamba markdupes is fragment-based and
# therefore superior to MACS2 mate 1-based duplicate detection


### MACS2 peak calling #########################################################

if pairedEnd:
    rule writeFragmentSize:
        input: "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv"
        output: "MACS2/fragmentSize.metrix.tsv"


    rule MACS2:
        input:
            chip = "filtered_bam/{chip_sample}.filtered.bam",
            control =
                lambda wildcards: "filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam" if get_control(wildcards.chip_sample)
                else [],
            insert_size_metrics = "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv"
        output:
            peaks = "MACS2/{chip_sample}.filtered.BAM_peaks.xls",
            peaksPE = "MACS2/{chip_sample}.filtered.BAMPE_peaks.xls"
        params:
            broad_calling =
                lambda wildcards: "--broad " if is_broad(wildcards.chip_sample) else "",
            control_param =
                lambda wildcards: " -c filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam" if get_control(wildcards.chip_sample)
                else "",
            genome_size = lambda wildcards: "-g "+genome_size if not cutntag else " ",
            ext_size =
                lambda wildcards: " --nomodel --extsize "+get_pe_frag_length(wildcards.chip_sample) if not cutntag else " ",
            peakCaller_options = lambda wildcards: peakCallerOptions if not cutntag else " -p 1e-5 ",
            bampe_options = lambda wildcards: BAMPEPeaks if not cutntag else " ",
            bam_options = lambda wildcards: BAMPeaks if not cutntag else " "
        log:
            out = "MACS2/logs/MACS2.{chip_sample}.filtered.out",
            err = "MACS2/logs/MACS2.{chip_sample}.filtered.err"
        benchmark:
            "MACS2/.benchmark/MACS2.{chip_sample}.filtered.benchmark"
        conda: CONDA_CHIPSEQ_ENV
        shell: """
            macs2 callpeak -t {input.chip} {params.control_param} \
                -f BAM \
                {params.bam_options}
                {params.genome_size} \
                {params.ext_size} \
                --keep-dup all \
                --outdir MACS2 \
                --name {wildcards.chip_sample}.filtered.BAM \
                {params.peakCaller_options} \
                {params.broad_calling} > {log.out} 2> {log.err}

            # also run MACS2 in paired-end mode BAMPE for comparison with single-end mode
            macs2 callpeak -t {input.chip} \
                {params.control_param} -f BAMPE \
                {params.bampe_options} \
                {params.peakCaller_options} \
                {params.genome_size} --keep-dup all \
                --outdir MACS2 --name {wildcards.chip_sample}.filtered.BAMPE \
                {params.broad_calling} > {log.out}.BAMPE 2> {log.err}.BAMPE
            """
else:
    rule MACS2:
        input:
            chip = "filtered_bam/{chip_sample}.filtered.bam",
            control =
                lambda wildcards: "filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam" if get_control(wildcards.chip_sample)
                else []
        output:
            peaks = "MACS2/{chip_sample}.filtered.BAM_peaks.xls",
        params:
            genome_size = int(genome_size),
            broad_calling =
                lambda wildcards: "--broad" if is_broad(wildcards.chip_sample)
                else "",
            control_param =
                lambda wildcards: "-c filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam" if get_control(wildcards.chip_sample)
                else "",
            frag_size=fragmentLength,
            mfold=mfold,
            qval_cutoff=qval
        log:
            out = "MACS2/logs/MACS2.{chip_sample}.filtered.out",
            err = "MACS2/logs/MACS2.{chip_sample}.filtered.err"
        benchmark:
            "MACS2/.benchmark/MACS2.{chip_sample}.filtered.benchmark"
        conda: CONDA_CHIPSEQ_ENV
        shell: """
            macs2 callpeak -t {input.chip} {params.control_param} -f BAM -g {params.genome_size} --qvalue {params.qval_cutoff} --keep-dup all --outdir MACS2 \
                --name {wildcards.chip_sample}.filtered.BAM --mfold {params.mfold} --extsize {params.frag_size}\
                {params.broad_calling} > {log.out} 2> {log.err}
            """



### MACS2 peak quality control #################################################

rule MACS2_peak_qc:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        xls = "MACS2/{sample}.filtered.BAM_peaks.xls"
    output:
        qc = "MACS2/{sample}.filtered.BAM_peaks.qc.txt"
    params:
        peaks =
            lambda wildcards: "MACS2/{}.filtered.BAM_peaks.broadPeak".format(wildcards.sample) if is_broad(wildcards.sample)
                              else "MACS2/{}.filtered.BAM_peaks.narrowPeak".format(wildcards.sample),
        genome_index = genome_index
    benchmark:
        "MACS2/.benchmark/MACS2_peak_qc.{sample}.filtered.benchmark"
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
        bam = "filtered_bam/{sample}.filtered.bam"
    output:
        bam = temp("filtered_bam/{sample}.namesorted.bam")
    log:
        "filtered_bam/logs/{sample}.namesort.err"
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
            bams=lambda wildcards: expand(os.path.join("filtered_bam", "{sample}.namesorted.bam"), sample=genrichDict[wildcards.group]),
            control = lambda wildcards: ["filtered_bam/"+get_control(x)+".namesorted.bam" for x in genrichDict[wildcards.group]] if chip_samples_w_ctrl else []
        output:
            "Genrich/{group}.narrowPeak"
        log: "Genrich/logs/{group}.log"
        params:
            bams = lambda wildcards: ",".join(expand(os.path.join("filtered_bam", "{sample}.namesorted.bam"), sample=genrichDict[wildcards.group])),
            blacklist = "-E {}".format(blacklist_bed) if blacklist_bed else "",
            control_pfx=lambda wildcards,input: "-c" if input.control else "",
            control=lambda wildcards,input: ",".join(input.control) if input.control else ""
        conda: CONDA_CHIPSEQ_ENV
        shell: """
            Genrich -t {params.bams} {params.control_pfx} {params.control} -o {output} -r {params.blacklist} -y 2> {log}
            """
else:
    rule Genrich_peaks:
        input:
            bams=lambda wildcards: expand(os.path.join("filtered_bam", "{sample}.namesorted.bam"), sample=genrichDict[wildcards.group]),
            control = lambda wildcards: ["filtered_bam/"+get_control(x)+".namesorted.bam" for x in genrichDict[wildcards.group] ] if chip_samples_w_ctrl else []
        output:
            "Genrich/{group}.narrowPeak"
        log: "Genrich/logs/{group}.log"
        params:
            bams = lambda wildcards: ",".join(expand(os.path.join("filtered_bam", "{sample}.namesorted.bam"), sample=genrichDict[wildcards.group])),
            blacklist = "-E {}".format(blacklist_bed) if blacklist_bed else "",
            control_pfx=lambda wildcards,input: "-c" if input.control else "",
            control=lambda wildcards,input: ",".join(input.control) if input.control else "",
            frag_size=fragmentLength
        conda: CONDA_CHIPSEQ_ENV
        shell: """
            Genrich -t {params.bams} {params.control_pfx} {params.control} -o {output} -r {params.blacklist} -w {params.frag_size} 2> {log}
            """
