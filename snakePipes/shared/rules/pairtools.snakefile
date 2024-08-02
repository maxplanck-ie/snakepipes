# based on https://github.com/caballero/snakemake-pairtools-phased/tree/df410ff

# Define function that returns pair files (phased or unphased), based on the reference.
def ret_pair(wildcards):
    if 'diploid_genome' in wildcards.ref:
        # Phased setting
        return f"pairs/{wildcards.sample}.{wildcards.ref}_phased.pairs.gz"
    else:
        return f"pairs/{wildcards.sample}.{wildcards.ref}.pairs.gz"

# different to bwa.snakefile
# here we skip the expensive sorting with samtools after bwa mem
# consider making this optional in bwa.snakefile
rule bwa_mapping:
    input:
        fq1 = "FASTQ_fastp/{sample}_R1.fastq.gz",
        fq2 = "FASTQ_fastp/{sample}_R2.fastq.gz",
        ix = "genome/{ref}.fa.gz.bwt"
    output:
        bam = "bam/{sample}.{ref}.bam"
    threads: 30
    params:
      bwathreads = config['alignerThreads'],
      bwaparams = config['alignerOptions'],
      fna = lambda wildcards, input: Path(input.ix).with_suffix('')
    resources:
      mem_mb = 3000
    conda: CONDA_MAKEPAIRS_ENV
    shell:
        """
        bwa mem \
            {params.bwaparams} \
            -t {params.bwathreads} \
            {params.fna} \
            {input.fq1} \
            {input.fq2} \
        | samtools view -@ 8 -b \
        > {output.bam}
        """

rule pairtools_parse:
    input:
        bam = "bam/{sample}.{ref}.bam",
        chr_sizes = "genome/{ref}.chromsizes"
    output:
        pairs = "pairs/{sample}.{ref}.pairs.gz"
    params:
        minmapq = 40,
        cols = lambda wildcards: '--add-columns XB,AS,XS' if 'diploid_genome' in wildcards.ref else ''
    threads: 12
    conda: CONDA_MAKEPAIRS_ENV
    shell:
        """
        pairtools parse \
            --min-mapq {params.minmapq} \
            {params.cols} \
            --drop-sam \
            --walks-policy 5unique \
            -c {input.chr_sizes} \
            {input.bam} \
            -o {output.pairs}
        """

rule pairtools_phase:
    input:
        pairs = "pairs/{sample}.diploid_genome.pairs.gz"
    output:
        pairs = "pairs/{sample}.diploid_genome_phased.pairs.gz"
    params:
        hap1 = strains[0],
        hap2 = strains[1]
    threads: 12
    conda: CONDA_MAKEPAIRS_ENV
    shell:
        """
        pairtools phase \
            --phase-suffixes _{params.hap1} _{params.hap2} \
            --tag-mode XB \
            --clean-output \
            {input.pairs} -o {output.pairs}
        """

rule pairtools_sort:
    input:
        ret_pair
    output:
        pairs = "pairs/{sample}.{ref}.pairs.sorted.gz"
    threads: 20
    conda: CONDA_MAKEPAIRS_ENV
    shell:
        """
        pairtools sort \
            {input} \
            -o {output.pairs} \
            --memory 20G
        """

rule pairtools_dedup:
    input:
        pairs = "pairs/{sample}.{ref}.pairs.sorted.gz"
    output:
        pairs = "pairs/{sample}.{ref}.pairs.dedup.gz",
        stats = "pairs/{sample}.{ref}.pairs.dedup.stats"
    params:
        extra_cols = lambda wildcards: '--extra-col-pair phase1 phase2'  if 'diploid' in wildcards.ref else ''
    threads: 12
    conda: CONDA_MAKEPAIRS_ENV
    shell:
        """
        pairtools dedup \
            --mark-dups \
            {params.extra_cols} \
            --output-dups - \
            --output-unmapped - \
            --output-stats {output.stats} \
            -o {output.pairs} \
            {input.pairs}
        """

rule pairtools_filter_phased:
    input:
        pairs = "pairs/{sample}.diploid_genome.pairs.dedup.gz"
    output:
        stats = "phase_stats/{sample}.diploid_genome_{phasetype}.pairs.stats",
        pairs = "phase_stats/{sample}.diploid_genome_{phasetype}.pairs.gz"
    params:
        filterparam = lambda wildcards: PHASEDIC[wildcards.phasetype]
    resources:
        mem_mb = 1000
    threads: 12
    conda: CONDA_MAKEPAIRS_ENV
    shell:
        """
        pairtools select \
            '{params.filterparam}' \
            {input.pairs} \
            -o {output.pairs}
        pairtools stats {output.pairs} -o {output.stats}
        """

rule multiqc:
    input:
        stats = expand(
            "pairs/{sample}.{ref}.pairs.dedup.stats",
            sample=samples,
            ref=REFERENCES
        ),
        phasedstats = expand(
            "phase_stats/{sample}.diploid_genome_{phasetype}.pairs.stats",
            sample=samples,
            phasetype = PHASEDIC.keys()
        )
    output:
        html = "multiqc/multiqc_report.html"
    params:
        odir = "multiqc"
    threads: 1
    conda: CONDA_MAKEPAIRS_ENV
    shell:
        """
        echo input: {input.phasedstats}
        multiqc \
            --module pairtools \
            --module fastqc \
            --module fastp \
            -o {params.odir} \
            .
        """
