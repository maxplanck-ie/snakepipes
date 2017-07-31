# rule filterMitochondrion:
#     input:
#         bam = "filtered_bam/{sample}.filtered.bam"
#     output:
#         bam = "filtered_bam/{sample}.filtered.noMitochondrion.bam"
#     shell:
#         # step1 identify all chroms and produce bed, w/o mitochondrion
#         samtools idxstats ATAC_S11.bam | cut -f 1 | grep 'mitochondrion'
#         # step2 write reads according to bed file

rule sortByName:
    input:
        "filtered_bam/{sample}.filtered.bam"
    output:
        "peaks_openChromatin/{sample}.filtered.sorted.bam"
    params:
        byQuery='-n'
    threads: 6
    log: "peaks_openChromatin/logs/sortByName/{sample}.log"
    shell:
        samtools_path+"samtools sort {params.byQuery} -@ {threads} {input} -o {output} &> {log}" ## TMPDIR (environment variable) for scratch usage

rule reads2fragments:
    input:
        "peaks_openChromatin/{sample}.filtered.sorted.bam"
    output:
        "peaks_openChromatin/{sample}.all.bedpe"
    shell:
        "/package/bedtools2/bin/bedtools bamtobed -bedpe -i {input} | "
        "awk -v OFS='\\t' -v pos_offset=\"4\" -v neg_offset=\"5\" "
        "'{{ print($1, $2 - pos_offset , $6 + neg_offset ) }}' > {output}"

rule filterFragments:
    input:
        "peaks_openChromatin/{sample}.all.bedpe"
    output:
        "peaks_openChromatin/{sample}.openchrom.bedpe"
    params:
        cutoff = atac_fragment_cutoff
    shell:
        "cat {input} | "
        "awk -v cutoff={params.cutoff} -v OFS='\\t' \"{{ if(\$3-\$2 < cutoff) {{ print (\$0) }} }}\" > "
        "{output}"

rule callOpenChromatin:
    input:
        "peaks_openChromatin/{sample}.openchrom.bedpe"
    output:
        peaks='peaks_openChromatin/openchromatin_{sample}_peaks.narrowPeak',
        pileup='peaks_openChromatin/openchromatin_{sample}_treat_pileup.bdg',
        ctrl='peaks_openChromatin/openchromatin_{sample}_control_lambda.bdg'
    params:
        directory = "peaks_openChromatin",
        genome=genome[0:2],
        name='openchromatin_{sample}',
        bandwidth='--bw 25', # + bw_binsize
        qval_cutoff='--qvalue 0.01',
        nomodel='--nomodel',
        write_bdg='--bdg',
        fileformat='--format BEDPE'
    threads: 6
    log: "peaks_openChromatin/logs/callOpenChromatin/{sample}_macs2.log"
    shell: # or run:
        ## macs2
        macs2_path+"macs2 callpeak "
            "--treatment {input} "
            "--gsize {params.genome} "
            "--name {params.name} "
            "--outdir {params.directory} "
            "{params.fileformat} {params.bandwidth} {params.qval_cutoff} {params.nomodel} {params.write_bdg} "
            "&> {log}"
