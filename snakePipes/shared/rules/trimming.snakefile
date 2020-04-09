### cutadapt #################################################################
if pairedEnd:
    rule cutadapt:
        input:
            r1 = fastq_indir_trim+"/{sample}"+reads[0]+".fastq.gz",
            r2 = fastq_indir_trim+"/{sample}"+reads[1]+".fastq.gz"
        output:
            r1 = "FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
            r2 = "FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz"
        params:
            opts = lambda wildcards: str(trimmerOptions or '')
        log:
            out = "FASTQ_Cutadapt/logs/Cutadapt.{sample}.out",
            err = "FASTQ_Cutadapt/logs/Cutadapt.{sample}.err"
        benchmark:
            "FASTQ_Cutadapt/.benchmark/Cutadapt.{sample}.benchmark"
        threads: 8
        conda: CONDA_SHARED_ENV
        shell: """
            cutadapt -j {threads} -e 0.1 -q 16 -O 3 --trim-n --minimum-length 25 -a AGATCGGAAGAGC -A AGATCGGAAGAGC {params.opts} \
                -o "{output.r1}" -p "{output.r2}" "{input.r1}" "{input.r2}" > {log.out} 2> {log.err}
            """
else:
    rule cutadapt:
        input:
            r1 = fastq_indir_trim+"/{sample}"+reads[0]+".fastq.gz",
        output:
            "FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
        params:
            opts = lambda wildcards: str(trimmerOptions or '')
        log:
            out = "FASTQ_Cutadapt/logs/Cutadapt.{sample}.out",
            err = "FASTQ_Cutadapt/logs/Cutadapt.{sample}.err"
        benchmark:
            "FASTQ_Cutadapt/.benchmark/Cutadapt.{sample}.benchmark"
        threads: 8
        conda: CONDA_SHARED_ENV
        shell: """
            cutadapt -j {threads} -e 0.1 -q 16 -O 3 --trim-n --minimum-length 25 -a AGATCGGAAGAGC {params.opts} \
                -o "{output}" "{input.r1}" > {log.out} 2> {log.err}
            """


### fastp #################################################################
# TODO: (1) ensure that multiQC sees the json files (2) remove reads[0] from the json file for MultiQC rendering
if pairedEnd:
    rule fastp:
        input:
            fastq_indir_trim+"/{sample}"+reads[0]+".fastq.gz",
            fastq_indir_trim+"/{sample}"+reads[1]+".fastq.gz"
        output:
            "FASTQ_fastp/{sample}"+reads[0]+".fastq.gz",
            "FASTQ_fastp/{sample}"+reads[1]+".fastq.gz",
            "FASTQ_fastp/{sample}fastp.json",
            "FASTQ_fastp/{sample}fastp.html"
        params:
            opts = lambda wildcards: str(trimmerOptions or '')
        log:
            out = "FASTQ_fastp/logs/fastp.{sample}.out",
            err = "FASTQ_fastp/logs/fastp.{sample}.err"
        benchmark:
            "FASTQ_fastp/.benchmark/fastp.{sample}.benchmark"
        threads: 8
        conda: CONDA_SHARED_ENV
        shell: """
            fastp -w {threads} -i "{input[0]}" -I "{input[1]}" -o "{output[0]}" -O "{output[1]}" -j "{output[2]}" -h "{output[3]}" {params.opts} > {log.out} 2> {log.err}
            """
else:
    rule fastp:
        input:
            fastq_indir_trim+"/{sample}"+reads[0]+".fastq.gz"
        output:
            "FASTQ_fastp/{sample}"+reads[0]+".fastq.gz",
            "FASTQ_fastp/{sample}fastp.json",
            "FASTQ_fastp/{sample}fastp.html"
        params:
            opts = lambda wildcards: str(trimmerOptions or '')
        log:
            out = "FASTQ_fastp/logs/fastp.{sample}.out",
            err = "FASTQ_fastp/logs/fastp.{sample}.err"
        benchmark:
            "FASTQ_fastp/.benchmark/fastp.{sample}.benchmark"
        threads: 8
        conda: CONDA_SHARED_ENV
        shell: """
            fastp -w {threads} -i "{input[0]}" -o "{output[0]}" -j "{output[1]}" -h "{output[2]}" {params.opts} > {log.out} 2> {log.err}
            """


### TrimGalore #################################################################

if pairedEnd:
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
            opts = lambda wildcards: str(trimmerOptions or '')
        log:
            out = "FASTQ_TrimGalore/logs/TrimGalore.{sample}.out",
            err = "FASTQ_TrimGalore/logs/TrimGalore.{sample}.err"
        benchmark:
            "FASTQ_TrimGalore/.benchmark/TrimGalore.{sample}.benchmark"
        conda: CONDA_SHARED_ENV
        shell: """
            trim_galore --output_dir FASTQ_TrimGalore --paired --stringency 3 {params.opts} "{input.r1}" "{input.r2}" > {log.out} 2> {log.err}
            mv "{params.tmp1}" "{output.r1}"
            mv "{params.tmp2}" "{output.r2}"
            """
else:
    rule TrimGalore:
        input:
            fastq_indir_trim+"/{sample}"+reads[0]+".fastq.gz"
        output:
            "FASTQ_TrimGalore/{sample}"+reads[0]+".fastq.gz"
        params:
            tmp = "FASTQ_TrimGalore/{sample}_trimmed.fq.gz",
            opts = lambda wildcards: str(trimmerOptions or '')
        log:
            out = "FASTQ_TrimGalore/logs/TrimGalore.{sample}.out",
            err = "FASTQ_TrimGalore/logs/TrimGalore.{sample}.err"
        benchmark:
            "FASTQ_TrimGalore/.benchmark/TrimGalore.{sample}.benchmark"
        conda: CONDA_SHARED_ENV
        shell: """
            trim_galore --output_dir FASTQ_TrimGalore --stringency 3 {params.opts} "{input}" > {log.out} 2> {log.err}
            mv "{params.tmp}" "{output}"
            """


### FastQC_on_trimmed #######################################################

if pairedEnd:
    rule FastQC_on_trimmed:
        input:
            fastq_dir+"/{sample}{read}.fastq.gz"
        output:
            "FastQC_trimmed/{sample}{read}_fastqc.html"
        log:
            out = "FastQC_trimmed/logs/FastQC_trimmed.{sample}{read}.out",
            err = "FastQC_trimmed/logs/FastQC_trimmed.{sample}{read}.err"
        benchmark:
            "FastQC_trimmed/.benchmark/FastQC_trimmed.{sample}{read}.benchmark"
        threads: 2
        conda: CONDA_SHARED_ENV
        shell: """
            fastqc -o FastQC_trimmed "{input}" > {log.out} 2> {log.err}
            """
else:
    rule FastQC_on_trimmed_SE:
        input:
            fastq_dir+"/{sample}"+reads[0]+".fastq.gz"
        output:
            "FastQC_trimmed/{sample}"+reads[0]+"_fastqc.html"
        log:
            out = "FastQC_trimmed/logs/FastQC_trimmed.{sample}"+reads[0]+".out",
            err = "FastQC_trimmed/logs/FastQC_trimmed.{sample}"+reads[0]+".err"
        benchmark:
            "FastQC_trimmed/.benchmark/FastQC_trimmed.{sample}"+reads[0]+".benchmark"
        threads: 2
        conda: CONDA_SHARED_ENV
        shell: """
            fastqc -o FastQC_trimmed "{input}" > {log.out} 2> {log.err}
            """
