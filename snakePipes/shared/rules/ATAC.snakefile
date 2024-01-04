#rule filterFragments:
#    input:
#        "filtered_bam/{sample}.filtered.bam"
#    output:
#        shortBAM = temp(os.path.join(short_bams, "{sample}.short.bam")),
#        metrics = os.path.join(short_bams, "{sample}.short.metrics")
#    log: os.path.join(short_bams, "logs/{sample}.filterFragments.log")
#    params:
#        maxFragmentSize=maxFragmentSize,
#        minFragmentSize=minFragmentSize
#    threads: 6
#    conda: CONDA_SHARED_ENV
#    shell: """
#        alignmentSieve --bam {input} \
#        --outFile {output.shortBAM} -p {threads} \
#        --filterMetrics {output.metrics} \
#        --maxFragmentLength {params.maxFragmentSize} \
#        --minFragmentLength {params.minFragmentSize} \
#        2> {log}
#        """


#rule filterMetricsToHtml:
#    input:
#        expand(os.path.join(short_bams, "{sample}.short.metrics"), sample=samples)
#    output:
#        QCrep='Filtering_metrics/Filtering_report.html'
#    log:
#        err="Filtering_metrics/logs/produce_report.err",
#        out="Filtering_metrics/logs/produce_report.out"
#    conda: CONDA_RMD_ENV
#    threads: 1
#    script: "../rscripts/ATACseq_QC_report_template.Rmd"

rule filterFragments:
    input:
        "filtered_bam/{sample}.filtered.bam"
    output:
        shortBAM = temp(os.path.join(short_bams, "{sample}.short.bam"))
    log: os.path.join(short_bams, "logs/{sample}.filterFragments.log")
    params:
        maxFragmentSize=maxFragmentSize,
        minFragmentSize=minFragmentSize
    threads: 6
    conda: CONDA_SAMBAMBA_ENV
    shell: """
        sambamba view -f bam -F "template_length >= {params.minFragmentSize} and template_length <= {params.maxFragmentSize} or template_length >= -{params.maxFragmentSize} and template_length <= {params.minFragmentSize}" -t {threads} -o {output.shortBAM} {input}
        2> {log}
        """



# MACS2 BAMPE fails if there is just one fragment mapped
rule filterCoveragePerScaffolds:
    input:
        bam = os.path.join(short_bams, "{sample}.short.bam")
    output:
        whitelist = os.path.join(short_bams, "{sample}.chrom.whitelist"),
        shortbai = temp(os.path.join(short_bams, "{sample}.short.bam.bai")),
        bam = os.path.join(short_bams, "{sample}.short.cleaned.bam"),
        bai = os.path.join(short_bams, "{sample}.short.cleaned.bam.bai")
    log: os.path.join(short_bams, "logs/{sample}.filterCoveragePerScaffolds.log")
    params:
        count_cutoff = int(fragmentCountThreshold) * 2 # must contain more than 2 reads, i.e. 1 fragment
    threads: 6
    conda: CONDA_SHARED_ENV
    shell: """
        samtools index -@ {threads} {input.bam} 2> {log}
        samtools idxstats {input.bam} | awk -v cutoff={params.count_cutoff} \'$3 > cutoff\' | cut -f 1 > {output.whitelist} 2>> {log}
        samtools view -@ {threads} -bo {output.bam} {input.bam} $(cat {output.whitelist} | paste -sd\' \') 2>> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}

        """


# MACS2 BAMPE filter: samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048
rule callOpenChromatin:
    input:
        os.path.join(short_bams, "{sample}.short.cleaned.bam")
    output:
        peaks = os.path.join(outdir_MACS2, '{sample}.filtered.short.BAM_peaks.narrowPeak'),
        xls = os.path.join(outdir_MACS2, '{sample}.filtered.short.BAM_peaks.xls')
    params:
        directory = outdir_MACS2,
        genome_size = int(genome_size),
        name='{sample}',
        qval_cutoff=qval,
        nomodel='--nomodel',
        write_bdg='--bdg',
        fileformat='--format BAMPE'
    threads: 6
    log:
        out = os.path.join(outdir_MACS2, "logs", "callOpenChromatin", "{sample}_macs2.out"),
        err = os.path.join(outdir_MACS2, "logs", "callOpenChromatin", "{sample}_macs2.err")
    conda: CONDA_ATAC_ENV
    shell: """
        macs2 callpeak --treatment {input} \
            -g {params.genome_size} \
            --name {params.name}.filtered.short.BAM \
            --outdir {params.directory} \
            {params.fileformat} \
            --qvalue {params.qval_cutoff} \
            {params.nomodel} \
            {params.write_bdg} > {log.out} 2> {log.err}
        """


rule tempChromSizes:
    input: genome_index
    output: temp("HMMRATAC/chrom_sizes")
    log: "HMMRATAC/logs/tempChromSizes.log"
    shell: """
        cut -f 1,2 {input} > {output} 2> {log}
        """


# TODO: -q MINMAPQ -Xmx value is currently hard-coded
# Actually uses 2-4 cores, even though there's no option for it!
# Requires PE data
rule HMMRATAC_peaks:
    input:
        "filtered_bam/{sample}.filtered.bam",
        "filtered_bam/{sample}.filtered.bam.bai",
        "HMMRATAC/chrom_sizes"
    output:
        "HMMRATAC/{sample}.log",
        "HMMRATAC/{sample}.model",
        "HMMRATAC/{sample}_peaks.gappedPeak",
        "HMMRATAC/{sample}_summits.bed",
        "HMMRATAC/{sample}_training.bed"
    log: "HMMRATAC/logs/{sample}.HMMRATAC_peaks.log"
    params:
        blacklist = "-e {}".format(blacklist_bed) if blacklist_bed else ""
    conda: CONDA_ATAC_ENV
    threads: 4
    shell: """
        HMMRATAC -Xmx10G -b {input[0]} -i {input[1]} -g {input[2]} {params.blacklist} -o HMMRATAC/{wildcards.sample} 2> {log}
        """

#Genrich requires namesorted bams
rule namesort_bams:
    input:
        bam = short_bams + "{sample}.short.cleaned.bam"
    output:
        bam = temp(short_bams + "{sample}.short.namesorted.bam")
    log:
        short_bams + "logs/{sample}.namesort.err"
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
if not isMultipleComparison:
    rule Genrich_peaks:
        input:
            bams=lambda wildcards: expand(short_bams + "{sample}.short.namesorted.bam", sample=genrichDict[wildcards.group])
        output:
            "Genrich/{group}.narrowPeak"
        log: "Genrich/logs/{group}.Genrich_peaks.log"
        params:
            bams = lambda wildcards: ",".join(expand(short_bams + "{sample}.short.namesorted.bam", sample=genrichDict[wildcards.group])),
            blacklist = "-E {}".format(blacklist_bed) if blacklist_bed else ""
        conda: CONDA_ATAC_ENV
        shell: """
            Genrich  -t {params.bams} -o {output} -r {params.blacklist} -j -y 2> {log}
            """

else:
    rule Genrich_peaks:
        input:
            bams=lambda wildcards: expand(short_bams + "{sample}.short.namesorted.bam", sample=genrichDict[wildcards.compGroup][wildcards.group]),
        output:
            "Genrich/{group}.{compGroup}.narrowPeak"
        log: "Genrich/logs/{group}.{compGroup}.log"
        params:
            bams = lambda wildcards: ",".join(expand(os.path.join(short_bams, "{sample}.short.namesorted.bam"), sample=genrichDict[wildcards.compGroup][wildcards.group])),
            blacklist = "-E {}".format(blacklist_bed) if blacklist_bed else "",
        conda: CONDA_ATAC_ENV
        shell: """
            Genrich -t {params.bams} -o {output} -r {params.blacklist} -j -y 2> {log}
            """
