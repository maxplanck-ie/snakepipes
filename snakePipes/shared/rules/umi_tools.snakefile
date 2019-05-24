##umi_tools###############
if umibarcode:
    if paired:
        rule umi_extract:
            input:
                r1 = "FASTQ/downsample_{sample}"+reads[0]+".fastq.gz" if downsample else "originalFASTQ/{sample}"+reads[0]+".fastq.gz",
                r2 = "FASTQ/downsample_{sample}"+reads[1]+".fastq.gz" if downsample else "originalFASTQ/{sample}"+reads[1]+".fastq.gz"
            output:
                r1 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
                r2 = "FASTQ/{sample}"+reads[1]+".fastq.gz"
            params:
                bcpattern = str(bcpattern)
            log:
                out = "FASTQ/logs/{sample}_log.out",
                err = "FASTQ/logs/{sample}_log.err"
            conda: CONDA_DNA_MAPPING_ENV
            shell:"""
                umi_tools extract -I {input.r1} --read2-in={input.r2} \
                --bc-pattern={params.bcpattern} --stdout={output.r1} \
                --read2-out={output.r2} -L {log.out} -E {log.err}
                """

    else:
        rule umi_extract:
            input:
                r1 = "FASTQ/downsample_{sample}"+reads[0]+".fastq.gz" if downsample else "originalFASTQ/{sample}"+reads[0]+".fastq.gz",
            output:
                r1 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
            params:
                bcpattern = str(bcpattern)
            log:
                out = "FASTQ/logs/{sample}_log.out",
                err = "FASTQ/logs/{sample}_log.err"
            conda:
                CONDA_DNA_MAPPING_ENV
            shell: """
                umi_tools extract -I {input.r1} --stdout={output.r1} \
                --bc-pattern={params.bcpattern} -L {log.out} -E {log.err}
                """

else:
    rule FASTQ1:
          input:
              "FASTQ/downsample_{sample}"+reads[0]+".fastq.gz" if downsample else "originalFASTQ/{sample}"+reads[0]+".fastq.gz"
          output:
              "FASTQ/{sample}"+reads[0]+".fastq.gz"
          shell:
              "( [ -f {output} ] || ln -s -r {input} {output} )"

    rule FASTQ2:
          input:
              "FASTQ/downsample_{sample}"+reads[1]+".fastq.gz" if downsample else "originalFASTQ/{sample}"+reads[1]+".fastq.gz"
          output:
              "FASTQ/{sample}"+reads[1]+".fastq.gz"
          shell:
              "( [ -f {output} ] || ln -s -r {input} {output} )"


if umidedup:
    rule filter_reads:
        input:
            bamfile = "filtered_bam/{sample}.filtered.tmp.bam",
            indexfile = "filtered_bam/{sample}.filtered.tmp.bam.bai"
        output:
            bamfile = "filtered_bam/{sample}.filtered.bam"
        params:
            umitools_options = str(umidedup_opts),
            umi_sep = str(umidedup_sep),
        log:
            out = "filtered_bam/logs/umi_dedup.{sample}.out",
            err = "filtered_bam/logs/umi_dedup.{sample}.err"
        conda:
            CONDA_DNA_MAPPING_ENV
        shell: """
            umi_tools dedup -I {input.bamfile} \
            -S {output.bamfile} -L {log.out} -E {log.err} \
            --umi-separator {params.umi_sep} \
            {params.umitools_options}
            """
else:
    rule filter_reads:
        input:
            bamfile = "filtered_bam/{sample}.filtered.tmp.bam",
        output:
            bamfile = "filtered_bam/{sample}.filtered.bam"
        conda:
            CONDA_DNA_MAPPING_ENV
        shell: """
            ln -s -r {input.bamfile} {output.bamfile}
            """

rule samtools_index_filtered:
    input:
        "filtered_bam/{sample}.filtered.bam"
    output:
        "filtered_bam/{sample}.filtered.bam.bai"
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"
