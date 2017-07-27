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
        bam = "filtered_bam/{sample}.filtered.bam"
    output:
        bam = "peaks_openChromatin/{sample}.filtered.sorted.bam"
    params:
        byQuery='-n'
    threads: 6
    log: "peaks_openChromatin/logs/sortByName/{sample}.log"
    shell:
        samtools_path+"samtools sort {params.byQuery} -@ {threads} {input} -o {output} &> {log}" ## TMPDIR (environment variable) for scratch usage


## actually wrong fragment shifting
rule reads2fragments:
    input:
        bam = "peaks_openChromatin/{sample}.filtered.sorted.bam"
    output:
        bedpe = "peaks_openChromatin/{sample}.all.bedpe"
    shell:
            "/package/bedtools2/bin/bedtools bamtobed -bedpe -i {input} | cut -f 1,2,6 | "
            "awk -v OFS='\\t' -v pos_offset=\"4\" -v neg_offset=\"5\" "
            "'{{ if( $9 == \"+\" ) {{ print($1, $2 + pos_offset , $3 - neg_offset ) }} "
            "else {{print($1, $2 - neg_offset , $3 + pos_offset) }} }}' > {output}"

rule filterFragments:
    input:
        bedpe="peaks_openChromatin/{sample}.all.bedpe"
    output:
        bedpe="peaks_openChromatin/{sample}.openchrom.bedpe"
    params:
        cutoff = atac_fragment_cutoff
    shell:
        "cat {input} | "
        "awk -v cutoff={params.cutoff} -v OFS='\\t' \"{{ if(\$3-\$2 < cutoff) {{ print (\$0) }} }}\" > "
        "{output}"

rule fragmentSizeDistribution:
    input:
        "peaks_openChromatin/{sample}.all.bedpe"
    output:
        "peaks_openChromatin/{sample}.all.fragdistr"
    shell:
        "cat {input} | awk '{{ print($3 - $2) }}' | sort -h | uniq -c > {output}"

rule callOpenChromatin:
    input:
        bedpe= "peaks_openChromatin/{sample}.openchrom.bedpe"
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

rule bedGraphToBigWig:
    input:
        "peaks_openChromatin/openchromatin_{sample}_{condition}.bdg"
    output:
        "peaks_openChromatin/openchromatin_{sample}_{condition}.bw"
    params:
        chromsize=genome_index
    log: "peaks_openChromatin/logs/bedGraphToBigWig/openChromatin_{sample}_{condition}.log"
    shell:
        "/package/UCSCtools/bedGraphToBigWig {input} {params.chromsize} {output} &> {log}"
