

### cutadapt #################################################################
if paired:
    rule cutadapt:
        input:
            r1 = fastq_indir_trim+"/{sample}"+reads[0]+".fastq.gz",
            r2 = fastq_indir_trim+"/{sample}"+reads[1]+".fastq.gz"
        output:
            r1 = "FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
            r2 = "FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz"
        params:
            opts = str(trim_options or '')
        log:
            "FASTQ_Cutadapt/logs/Cutadapt.{sample}.log"
        benchmark:
            "FASTQ_Cutadapt/.benchmark/Cutadapt.{sample}.benchmark"
        threads: 8
        conda: CONDA_SHARED_ENV
        shell: """
            cutadapt {params.opts} -j {threads} -e 0.1 -q 16 -O 3 --trim-n --minimum-length 25 -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
                -o {output.r1} -p {output.r2} {input.r1} {input.r2} &> {log}
            """
else:
    rule cutadapt:
        input:
            r1 = fastq_indir_trim+"/{sample}.fastq.gz",
        output:
            "FASTQ_Cutadapt/{sample}.fastq.gz",
        params:
            opts = str(trim_options or '')
        log:
            "FASTQ_Cutadapt/logs/Cutadapt.{sample}.log"
        benchmark:
            "FASTQ_Cutadapt/.benchmark/Cutadapt.{sample}.benchmark"
        threads: 8
        conda: CONDA_SHARED_ENV
        shell: """
            cutadapt {params.opts} -j {threads} -e 0.1 -q 16 -O 3 --trim-n --minimum-length 25 -a AGATCGGAAGAGC \
                -o {output} {input.r1} &> {log}
            """


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
            opts = str(trim_options or '')
        log:
            "FASTQ_TrimGalore/logs/TrimGalore.{sample}.log"
        benchmark:
            "FASTQ_TrimGalore/.benchmark/TrimGalore.{sample}.benchmark"
        conda: CONDA_SHARED_ENV
        shell: """
            trim_galore --output_dir FASTQ_TrimGalore --paired --stringency 3 {params.opts} {input.r1} {input.r2} &> {log} && \
            (mv {params.tmp1} {output.r1} ; mv {params.tmp2} {output.r2})
            """
else:
    rule TrimGalore:
        input:
            fastq_indir_trim+"/{sample}.fastq.gz"
        output:
            "FASTQ_TrimGalore/{sample}.fastq.gz"
        params:
            tmp = "FASTQ_TrimGalore/{sample}_trimmed.fq.gz",
            opts = str(trim_options or '')
        log:
            "FASTQ_TrimGalore/logs/TrimGalore.{sample}.log"
        benchmark:
            "FASTQ_TrimGalore/.benchmark/TrimGalore.{sample}.benchmark"
        conda: CONDA_SHARED_ENV
        shell: """
            trim_galore --output_dir FASTQ_TrimGalore --stringency 3 {params.opts} {input} &> {log} && \
            mv {params.tmp} {output}
            """


### FastQC_on_trimmed #######################################################

rule FastQC_on_trimmed:
    input:
        fastq_dir+"/{sample}{read}.fastq.gz"
    output:
        "FastQC_trimmed/{sample}{read}_fastqc.html"
    log:
        "FastQC_trimmed/logs/FastQC_trimmed.{sample}{read}.log"
    benchmark:
        "FastQC_trimmed/.benchmark/FastQC_trimmed.{sample}{read}.benchmark"
    threads: 2
    conda: CONDA_SHARED_ENV
    shell: "fastqc -o FastQC_trimmed {input} &> {log}"
