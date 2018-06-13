rule filterFragments:
    input:
        "filtered_bam/{sample}.filtered.bam"
    output:
        shortBAM = temp(os.path.join(outdir_MACS2, "{sample}.short.bam")),
        metrics = os.path.join(outdir_MACS2, "{sample}.short.metrics")
    params:
        cutoff = atac_fragment_cutoff
    threads: 6
    conda: CONDA_SHARED_ENV
    shell: """
        alignmentSieve --bam {input} \
        --outFile {output.shortBAM} -p {threads} \
        --filterMetrics {output.metrics} \
        --maxFragmentLength {params.cutoff}
        """

# MACS2 BAMPE filter: samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048
rule callOpenChromatin:
    input:
        os.path.join(outdir_MACS2, "{sample}.short.bam")
    output:
        peaks = os.path.join(outdir_MACS2, '{sample}.filtered.BAM_peaks.narrowPeak'),
        xls = os.path.join(outdir_MACS2, '{sample}.filtered.BAM_peaks.xls')
    params:
        directory = outdir_MACS2,
        genome_size = int(genome_size),
        name='{sample}',
        qval_cutoff='--qvalue 0.001',
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
            --name {params.name}.filtered.BAM \
            --outdir {params.directory} \
            {params.fileformat} \
            {params.qval_cutoff} \
            {params.nomodel} \
            {params.write_bdg} > {log.out} 2> {log.err}
        """
