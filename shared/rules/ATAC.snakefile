rule reads2Frags:
    input:
        "filtered_bam/{sample}.filtered.bam"
    output:
        allFrags=os.path.join(outdir_MACS2, "{sample}.all.bedpe"),
        shortFrags=os.path.join(outdir_MACS2, "{sample}.short.bedpe")
    params:
        cutoff=atac_fragment_cutoff
    threads: 6
    shell:
        samtools_path + "samtools sort -l 0 -n -@ {threads} {input} | "
         + bedtools_path +"bedtools bamtobed -bedpe -i - |"
         "awk -v OFS='\\t' -v pos_offset=\"4\" -v neg_offset=\"5\" "
         "'{{ print($1, $2 - pos_offset , $6 + neg_offset ) }}' | tee \"{output.allFrags}\" |"
         "awk -v cutoff={params.cutoff} -v OFS='\\t' \"{{ if(\$3-\$2 < cutoff) {{ print (\$0) }} }}\" > "
         "{output.shortFrags}"

# MACS2 BAMPE filter: samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048
rule callOpenChromatin:
    input:
        rules.reads2Frags.output.shortFrags
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
        fileformat='--format BEDPE'
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
