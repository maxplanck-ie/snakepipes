### add barcodes from R1 to R2 #########

rule fastq_barcode:
        input:
            R1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
            R2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz"
        output:
            R2_barcoded = "FASTQ_barcoded/{sample}.barcoded.fastq.gz"
        params:
            UMI_length = 6,
            CELLI_length = 6
        benchmark:
            "fastq_barcoded/.benchmark/fastq_barcoded.{sample}.benchmark"
        threads: 2
        shell:"""
            paste <(paste - - - - < <(zcat {input.R1}))   <(paste - - - - < <(zcat {input.R2}))  | \
            tr '\t' '\n' | \
            awk -v CBAR_LEN={params.CELLI_length} -v UMI_LEN={params.UMI_length} \ 
            '{{print $0}}' | pigz -p 2 > {output.R2_barcoded}
            """
