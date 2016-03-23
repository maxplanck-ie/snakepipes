if paired:
    rule TrimGalore:
        input:
            r1 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
            r2 = "FASTQ/{sample}"+reads[1]+".fastq.gz"
        output:
            r1 = "TrimGalore/{sample}"+reads[0]+".fastq.gz",
            r2 = "TrimGalore/{sample}"+reads[1]+".fastq.gz"
        params:
            tmp1 = "TrimGalore/{sample}"+reads[0]+"_val_1.fq.gz",
            tmp2 = "TrimGalore/{sample}"+reads[1]+"_val_2.fq.gz",
            opts = trim_galore_opts
        log:    "TrimGalore/log/{sample}.TrimGalore.log"
        benchmark:  "TrimGalore/.benchmark/TrimGalore.{sample}.benchmark"
        shell:
            os.path.join(cutadapt_dir,'activate')+" && "
            ""+trim_galore_path+" --paired {params.opts} --output_dir TrimGalore {input.r1} {input.r2} "
            "&& (mv {params.tmp1} {output.r1}; mv {params.tmp2} {output.r2}) "
            "2>&1 | tee {log}"
else:
    rule TrimGalore:
        input:  "FASTQ/{sample}.fastq.gz"
        output: "TrimGalore/{sample}.fastq.gz"
        params:
            tmp = "TrimGalore/{sample}_trimmed.fq.gz",
            opts = trim_galore_opts
        log:    "TrimGalore/log/{sample}.TrimGalore.log"
        benchmark:  "TrimGalore/.benchmark/TrimGalore.{sample}.benchmark"
        shell:
            os.path.join(cutadapt_dir,'activate')+" && "
            ""+trim_galore_path+" {params.opts} --output_dir TrimGalore {input} "
            "&& mv {params.tmp} {output} "
            "2>&1 | tee {log}"


rule FastQC_on_TrimGalore:
    input:  "TrimGalore/{sample}{read}.fastq.gz"
    output: "FastQC_on_TrimGalore/{sample}{read}_fastqc.html"
    log:    "FastQC_on_TrimGalore/log/{sample}{read}.FastQC_on_TrimGalore.log"
    benchmark:  "FastQC_on_TrimGalore/.benchmark/FastQC_on_TrimGalore.{sample}{read}.benchmark"
    shell:  fastqc_path+" -o FastQC_on_TrimGalore {input} 2>&1 | tee {log}"
