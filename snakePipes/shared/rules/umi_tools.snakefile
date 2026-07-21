##umi_tools###############

if UMIBarcode:
    if pairedEnd:
        rule umi_extract:
            input:
                r1 = "originalFASTQ/downsample_{sample}"+reads[0]+".fastq.gz" if downsample else "originalFASTQ/{sample}"+reads[0]+".fastq.gz",
                r2 = "originalFASTQ/downsample_{sample}"+reads[1]+".fastq.gz" if downsample else "originalFASTQ/{sample}"+reads[1]+".fastq.gz"
            output:
                r1 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
                r2 = "FASTQ/{sample}"+reads[1]+".fastq.gz"
            params:
                bcpattern = str(bcPattern)
            conda: CONDA_SHARED_ENV
            shell:"""
                umi_tools extract -I {input.r1} --read2-in={input.r2} \
                --bc-pattern={params.bcpattern} --bc-pattern2={params.bcpattern}\
                --stdout={output.r1} \
                --read2-out={output.r2}
                """

    else:
        rule umi_extract:
            input:
                r1 = "originalFASTQ/downsample_{sample}"+reads[0]+".fastq.gz" if downsample else "originalFASTQ/{sample}"+reads[0]+".fastq.gz",
            output:
                r1 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
            params:
                bcpattern = str(bcPattern)
            conda: CONDA_SHARED_ENV
            shell: """
                umi_tools extract -I {input.r1} --stdout={output.r1} \
                --bc-pattern={params.bcpattern}
                """

else:
    rule FASTQ1:
          input:
              "originalFASTQ/downsample_{sample}"+reads[0]+".fastq.gz" if downsample else "originalFASTQ/{sample}"+reads[0]+".fastq.gz"
          output:
              "FASTQ/{sample}"+reads[0]+".fastq.gz"
          shell: """
                ln -s ../{input} {output}
          """

    if pairedEnd or pipeline=="scrnaseq":
        rule FASTQ2:
            input:
                "originalFASTQ/downsample_{sample}"+reads[1]+".fastq.gz" if downsample else "originalFASTQ/{sample}"+reads[1]+".fastq.gz"
            output:
                "FASTQ/{sample}"+reads[1]+".fastq.gz"
            shell: """
                ln -s ../{input} {output}
          """

#If DNAmapping:
if UMIDedup:
    rule filter_reads_umi:
        input:
            bamfile = "filtered_bam/{sample}.filtered.tmp.bam",
            indexfile = "filtered_bam/{sample}.filtered.tmp.bam.bai"
        output:
            bamfile = "filtered_bam/{sample}.filtered.bam"
        params:
            umitools_options = str(UMIDedupOpts or ''),
            umitools_paired = "--paired " if pairedEnd else " ",
            umi_sep = str(UMIDedupSep),
        conda: CONDA_SHARED_ENV
        shell: """
            umi_tools dedup -I {input.bamfile} \
            -S {output.bamfile} \
            --umi-separator {params.umi_sep} \
            {params.umitools_paired} {params.umitools_options}
            """
else:
    rule filter_reads:
        input:
            bamfile = "filtered_bam/{sample}.filtered.tmp.bam"
        output:
            bamfile = "filtered_bam/{sample}.filtered.bam"
        shell: """
            mv -v {input} {output}
          """

rule samtools_index_filtered:
    input:
        "filtered_bam/{sample}.filtered.bam"
    output:
        "filtered_bam/{sample}.filtered.bam.bai"
    conda: CONDA_SHARED_ENV
    shell: "samtools index {input}"
