### TrimGalore #################################################################

if paired:
    rule TrimGalore:
        input:
            r1 = fastq_indir_trim+"/{sample}"+reads[0]+".fastq.gz",
            r2 = fastq_indir_trim+"/{sample}"+reads[1]+".fastq.gz"
        output:
            r1 = "FASTQ_TrimGalore/{sample}"+reads[0]+".fastq.gz",
            r2 = "FASTQ_TrimGalore/{sample}"+reads[1]+".fastq.gz"
        params:
            tmp1 = "FASTQ_TrimGalore/{sample}"+reads[0]+"_val_1.fq.gz",
            tmp2 = "FASTQ_TrimGalore/{sample}"+reads[1]+"_val_2.fq.gz",
            opts = trim_galore_opts
        log:
            "FASTQ_TrimGalore/logs/TrimGalore.{sample}.log"
        benchmark:
            "FASTQ_TrimGalore/.benchmark/TrimGalore.{sample}.benchmark"
        shell:
            trim_galore_path+"trim_galore "
                "--path_to_cutadapt "+cutadapt_path+"cutadapt "
                "--output_dir FASTQ_TrimGalore "
                " --paired "
                "{params.opts} {input.r1} {input.r2} "
                "&> {log} "
                "&& (mv {params.tmp1} {output.r1}; mv {params.tmp2} {output.r2})"
else:
    rule TrimGalore:
        input:
            fastq_indir_trim+"/{sample}.fastq.gz"
        output:
            "FASTQ_TrimGalore/{sample}.fastq.gz"
        params:
            tmp = "FASTQ_TrimGalore/{sample}_trimmed.fq.gz",
            opts = trim_galore_opts
        log:
            "FASTQ_TrimGalore/logs/TrimGalore.{sample}.log"
        benchmark:
            "FASTQ_TrimGalore/.benchmark/TrimGalore.{sample}.benchmark"
        shell:
            trim_galore_path+"trim_galore "
                "--path_to_cutadapt "+cutadapt_path+"cutadapt "
                "--output_dir FASTQ_TrimGalore "
                "{params.opts} {input} "
                "&> {log} "
                "&& mv {params.tmp} {output}"


### FastQC_on_TrimGalore #######################################################

rule FastQC_on_TrimGalore:
    input:
        "FASTQ_TrimGalore/{sample}{read}.fastq.gz"
    output:
        "FastQC_on_TrimGalore/{sample}{read}_fastqc.html"
    log:
        "FastQC_on_TrimGalore/logs/FastQC_on_TrimGalore.{sample}{read}.log"
    benchmark:
        "FastQC_on_TrimGalore/.benchmark/FastQC_on_TrimGalore.{sample}{read}.benchmark"
    threads: 2
    shell:
        fastqc_path+"fastqc -o FastQC_on_TrimGalore {input} &> {log}"
