rule generate_chromsizes:
  input:
    genome_index
  output:
    'genome/genome.chromsizes'
  threads: 1
  shell:'''
  cut -f1,2 {input} > {output}
  '''

rule bwa_mapping:
    input:
        fq1="FASTQ_fastp/{sample}_R1.fastq.gz",
        fq2="FASTQ_fastp/{sample}_R2.fastq.gz"
    output:
        bam="bam/{sample}.bam",
    threads: 30
    params:
        bwaparams=config["alignerOptions"],
        bwa_index = bwa_index
    resources:
        mem_mb=3000,
    benchmark:
        "bam/.benchmark/bwa_mapping.{sample}.benchmark"
    conda:
        CONDA_MAKEPAIRS_ENV
    shell:
        """
        bwa mem \
            {params.bwaparams} \
            -t 22 \
            {params.bwa_index} \
            {input.fq1} \
            {input.fq2} \
        | samtools view -@ 8 -b \
        > {output.bam}
        """

rule pairtools_parse:
    input:
        bam="bam/{sample}.bam",
        chr_sizes='genome/genome.chromsizes'
    output:
        pairs=temp("pairs/{sample}.unsorted.pairs.gz"),
    params:
        minmapq=40
    threads: 12
    benchmark:
        "pairs/.benchmark/pairtools_parse.{sample}.benchmark"
    conda:
        CONDA_MAKEPAIRS_ENV
    shell:
        """
        pairtools parse \
            --min-mapq {params.minmapq} \
            --drop-sam \
            --walks-policy 5unique \
            -c {input.chr_sizes} \
            {input.bam} \
            -o {output.pairs}
        """

rule pairtools_sort:
    input:
        pairs = "pairs/{sample}.unsorted.pairs.gz",
    output:
        pairs = "pairs/{sample}.pairs.gz",
    threads: 20
    benchmark:
        "pairs/.benchmark/pairtools_sort.{sample}.benchmark"
    conda:
        CONDA_MAKEPAIRS_ENV
    shell:
        """
        pairtools sort \
            {input.pairs} \
            -o {output.pairs} \
            --memory 20G
        """

rule pairtools_dedup:
    input:
        pairs="pairs/{sample}.pairs.gz",
    output:
        pairs="pairs/{sample}.pairs.dedup.gz",
        stats="pairs/{sample}.pairs.dedup.stats"
    threads: 12
    benchmark:
        "pairs/.benchmark/pairtools_dedup.{sample}.benchmark"
    conda:
        CONDA_MAKEPAIRS_ENV
    shell:
        """
        pairtools dedup \
            --mark-dups \
            --output-dups - \
            --output-unmapped - \
            --output-stats {output.stats} \
            -o {output.pairs} \
            {input.pairs}
        """

rule pairix:
    input:
        pairs = 'pairs/{sample}.pairs.dedup.gz'
    output:
        ix = 'pairs/{sample}.pairs.dedup.gz.px2'
    threads: 2
    conda:
        CONDA_MAKEPAIRS_ENV
    shell:
        """
        pairix -f -p pairs {input.pairs}
        """

rule cooler:
    input:
        pairs = 'pairs/{sample}.pairs.dedup.gz',
        ix = 'pairs/{sample}.pairs.dedup.gz.px2',
        chromsizes = 'genome/genome.chromsizes'
    output:
        cool = 'cooler/{sample}.5000.cool'
    threads: 20
    conda:
        CONDA_MAKEPAIRS_ENV
    shell:
        """
        cooler cload pairix -p {threads} {input.chromsizes}:5000 {input.pairs} {output.cool}
        cooler balance --nproc {threads} {output.cool}
        """

rule mcool:
    input:
        cool = 'cooler/{sample}.5000.cool'
    output:
        mcool = 'cooler/{sample}.5000.mcool'
    threads: 20
    conda:
        CONDA_MAKEPAIRS_ENV
    shell:
        """
        cooler zoomify --resolutions 5000,10000,20000,40000,80000,120000 --balance --nproc {threads} {input.cool}
        """

rule multiqc:
    input:
        cools=expand(
            'cooler/{sample}.5000.mcool', sample=samples
        )
    output:
        html="multiqc/multiqc_report.html",
    params:
        odir="multiqc",
    benchmark:
        "multiqc/.benchmark/multiqc.benchmark"
    threads: 1
    conda:
        CONDA_MAKEPAIRS_ENV
    shell:
        """
        multiqc \
            --module pairtools \
            --module fastqc \
            --module fastp \
            -o {params.odir} \
            .
        """