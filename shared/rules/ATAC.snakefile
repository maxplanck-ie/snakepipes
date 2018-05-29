rule filterFragments:
    input:
        "filtered_bam/{sample}.filtered.bam"
    output:
        shortBAM = os.path.join(outdir_MACS2, "{sample}.short.bam"),
        metrics = os.path.join(outdir_MACS2, "{sample}.short.metrics")
    params:
        cutoff = atac_fragment_cutoff
    threads: 6
    conda: CONDA_SHARED_ENV
    shell: """
        alignmentSieve --bam {input} --outFile {output.shortBAM} -p {threads} \
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
        bandwidth='--bw 25',
        qval_cutoff='--qvalue 0.01',
        nomodel='--nomodel',
        write_bdg='--bdg',
        fileformat='--format BAMPE'
    threads: 6
    log: os.path.join(outdir_MACS2, "logs", "callOpenChromatin","{sample}_macs2.log")
    conda: CONDA_SHARED_ENV
    shell: """
        macs2 callpeak --treatment {input} \
            -g {params.genome_size} \
            --name {params.name}.filtered.BAM \
            --outdir {params.directory} \
            --slocal 10000 \
            --nolambda \
            {params.fileformat} {params.bandwidth} {params.qval_cutoff} {params.nomodel} {params.write_bdg} &> {log}
        """
