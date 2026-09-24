### fastp #################################################################
# TODO: (1) ensure that multiQC sees the json files (2) remove reads[0] from the json file for MultiQC rendering
if pairedEnd:
    rule fastp:
        input:
            fastq_indir_trim+"/{sample}"+reads[0]+".fastq.gz",
            fastq_indir_trim+"/{sample}"+reads[1]+".fastq.gz"
        output:
            "FASTQ_fastp/{sample}"+reads[0].replace("R1", "")+".fastq.gz",
            "FASTQ_fastp/{sample}fastp.json",
            "FASTQ_fastp/{sample}fastp.html"
        params:
            opts = lambda wildcards: str(trimmerOptions or '')
        benchmark:
            "FASTQ_fastp/.benchmark/fastp.{sample}.benchmark"
        threads: lambda wildcards: 8 if 8<max_thread else max_thread
        conda: CONDA_SHARED_ENV
        shell: """
            fastp -w {threads} -i "{input[0]}" -I "{input[1]}" --adapter_fasta {Adapterseq} "--include_unmerged"  "{output[0]}" -j "{output[2]}" -h "{output[3]}" {params.opts}
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
        benchmark:
            "FASTQ_fastp/.benchmark/fastp.{sample}.benchmark"
        threads: lambda wildcards: 8 if 8<max_thread else max_thread
        conda: CONDA_SHARED_ENV
        shell: """
            fastp -w {threads} -i "{input[0]}" --adapter_fasta {Adapterseq} -o "{output[0]}" -j "{output[1]}" -h "{output[2]}" {params.opts}
            """

### FastQC_on_trimmed #######################################################

if pairedEnd:
    rule FastQC_on_trimmed:
        input:
            fastq_dir+"/{sample}{read}.fastq.gz"
        output:
            "FastQC_trimmed/{sample}{read}_fastqc.html"
        benchmark:
            "FastQC_trimmed/.benchmark/FastQC_trimmed.{sample}{read}.benchmark"
        threads: lambda wildcards: 2 if 2<max_thread else max_thread
        conda: CONDA_SHARED_ENV
        shell: """
            fastqc -o FastQC_trimmed "{input}"
            """
else:
    rule FastQC_on_trimmed_SE:
        input:
            fastq_dir+"/{sample}"+reads[0]+".fastq.gz"
        output:
            "FastQC_trimmed/{sample}"+reads[0]+"_fastqc.html"
        benchmark:
            "FastQC_trimmed/.benchmark/FastQC_trimmed.{sample}"+reads[0]+".benchmark"
        threads: lambda wildcards: 2 if 2<max_thread else max_thread
        conda: CONDA_SHARED_ENV
        shell: """
            fastqc -o FastQC_trimmed "{input}"
            """
