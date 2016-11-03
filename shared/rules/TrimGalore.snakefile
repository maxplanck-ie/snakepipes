### TrimGalore #################################################################
if paired:
    rule cutadapt:
        input:
            r1 = fastq_indir_trim+"/{sample}"+reads[0]+".fastq.gz",
            r2 = fastq_indir_trim+"/{sample}"+reads[1]+".fastq.gz"
        output:
            r1 = "FASTQ_Cutadapt/{sample}"+reads[0]+".fastq.gz",
            r2 = "FASTQ_Cutadapt/{sample}"+reads[1]+".fastq.gz"
        params:
            tmp1 = "{sample}"+reads[0]+".fq.gz",
            tmp2 = "{sample}"+reads[1]+".fq.gz",
            opts = trim_galore_opts
        log:
            "FASTQ_Cutadapt/logs/Cutadapt.{sample}.log"
        benchmark:
            "FASTQ_Cutadapt/.benchmark/Cutadapt.{sample}.benchmark"
        shell:
            cutadapt_path+"cutadapt "
                "-f fastq -e 0.1 -q 20 -O 2 --trim-n --minimum-length 25 -a AGATCGGAAGAGC -A AGATCGGAAGAGC "
                "-o ${{TMPDIR}}{params.tmp1} -p ${{TMPDIR}}{params.tmp2} " 
                "{input.r1} {input.r2} "
                "&> {log} "
                "&& (mv ${{TMPDIR}}{params.tmp1} {output.r1}; mv ${{TMPDIR}}{params.tmp2} {output.r2}; touch {output.r1} {output.r2})"
else:
    rule cutadapt:
        input:
            r1 = fastq_indir_trim+"/{sample}.fastq.gz",
        output:
            r1 = "FASTQ_Cutadapt/{sample}.fastq.gz",
        params:
            tmp = "{sample}.fq.gz",
            opts = trim_galore_opts
        log:
            "FASTQ_Cutadapt/logs/Cutadapt.{sample}.log"
        benchmark:
            "FASTQ_Cutadapt/.benchmark/Cutadapt.{sample}.benchmark"
        shell:
            cutadapt_path+"cutadapt "
                "-f fastq -e 0.1 -q 20 -O 2 --trim-n --minimum-length 25 -a AGATCGGAAGAGC "
                "-o ${{TMPDIR}}{params.tmp} " 
                "{input.r1} "
                "&> {log} "
                "&& (mv ${{TMPDIR}}{params.tmp} {output.r1}; touch {output.r1})"

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
    shell:
        fastqc_path+"fastqc -o FastQC_trimmed {input} &> {log}"
