
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
        os.path.join(outdir_MACS2, "{sample}.filtered.sorted.bam")
    params:
        byQuery='-n'
    threads: 6
    log: os.path.join(outdir_MACS2, "logs","sortByName","{sample}.log")
    shell:
        samtools_path+"samtools sort {params.byQuery} -@ {threads} {input} -o {output} &> {log}" ## TMPDIR (environment variable) for scratch usage

rule reads2fragments:
    input:
        os.path.join(outdir_MACS2, "{sample}.filtered.sorted.bam")
    output:
        os.path.join(outdir_MACS2, "{sample}.all.bedpe")
    params:
        samtools_macs2filter="-f 2 -F 4 -F 8 -F 256 -F 512 -F 2048"
    shell:
        samtools_path + "samtools view -b {params.samtools_macs2filter} {input} |"
        "/package/bedtools2/bin/bedtools bamtobed -bedpe -i - | "
        "awk -v OFS='\\t' -v pos_offset=\"4\" -v neg_offset=\"5\" "
        "'{{ print($1, $2 - pos_offset , $6 + neg_offset ) }}' > {output}"

rule filterChromosomes:
    input:
        os.path.join(outdir_MACS2, "{sample}.all.bedpe")
    output:
        os.path.join(outdir_MACS2, "{sample}.all_filtered.bedpe")
    params:
        chromlist = '^'+'(' + '|'.join( ('dmel_mitochondrion_genome', '3R','3L','2R', '2L', 'X', '4', 'Y') ) + ')' + '[[:space:]]' # dm6!
    shell:
        "egrep \'{params.chromlist}\' {input} > {output}"


rule filterNucleosomalFragments:
    input:
        os.path.join(outdir_MACS2, "{sample}.all_filtered.bedpe")
    output:
        os.path.join(outdir_MACS2, "{sample}.openchrom.bedpe")
    params:
        cutoff = atac_fragment_cutoff
    shell:
        "cat {input} | "
        "awk -v cutoff={params.cutoff} -v OFS='\\t' \"{{ if(\$3-\$2 < cutoff) {{ print (\$0) }} }}\" > "
        "{output}"


# samtools view -b -f 2 -F 4 -F 8 -F 256 -F 512 -F 2048

rule callOpenChromatin:
    input:
        os.path.join(outdir_MACS2, "{sample}.openchrom.bedpe")
    output:
        peaks = os.path.join(outdir_MACS2, 'openchromatin_{sample}_peaks.narrowPeak'),
        pileup = os.path.join(outdir_MACS2, 'openchromatin_{sample}_treat_pileup.bdg'),
        ctrl = os.path.join(outdir_MACS2, 'openchromatin_{sample}_control_lambda.bdg')
    params:
        directory = outdir_MACS2,
        genome=genome[0:2],
        name='openchromatin_{sample}',
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
            "{params.fileformat} {params.bandwidth} {params.qval_cutoff} {params.nomodel} {params.write_bdg} "
            "&> {log}"
