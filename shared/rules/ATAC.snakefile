rule filterFragments:
    input:
        "filtered_bam/{sample}.filtered.bam"
    output:
         shortFrags=os.path.join(outdir_MACS2, "{sample}.short.bam")
    params:
        cutoff=atac_fragment_cutoff,
        metrics=os.path.join(outdir_MACS2, "{sample}.short.metrics")
    threads: 6
    shell:
        deepTools_path + "alignmentSieve --bam {input} --outFile {output} "
        "--numberOfProcessors {threads} "
        "--filterMetrics {params.metrics} "
        "--maxFragmentLength {params.cutoff} "

# MACS2 BAMPE filter: samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048
rule callOpenChromatin:
    input:
        rules.filterFragments.output
    output:
        peaks = os.path.join(outdir_MACS2, '{sample}_peaks.narrowPeak'),
        pileup = os.path.join(outdir_MACS2, '{sample}_treat_pileup.bdg'),
        ctrl = os.path.join(outdir_MACS2, '{sample}_control_lambda.bdg'),
        xls = os.path.join(outdir_MACS2, '{sample}_peaks.xls')
    params:
        directory = outdir_MACS2,
        genome=genome[0:2],
        name='{sample}',
        bandwidth='--bw 25', # + bw_binsize
        qval_cutoff='--qvalue 0.01',
        nomodel='--nomodel',
        write_bdg='--bdg',
        fileformat='--format BAMPE'
    threads: 6
    log: os.path.join(outdir_MACS2, "logs", "callOpenChromatin","{sample}_macs2.log")
    shell: # or run:
        ## macs2
        macs2_path+"macs2 callpeak "
            "--treatment {input} "
            "--gsize {params.genome} "
            "--name {params.name} "
            "--outdir {params.directory} "
            "--slocal 10000 "
            "--nolambda "
            "{params.fileformat} {params.bandwidth} {params.qval_cutoff} {params.nomodel} {params.write_bdg} "
            "&> {log}"
