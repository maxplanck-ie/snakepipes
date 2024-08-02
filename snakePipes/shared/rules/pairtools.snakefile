# Define function that returns pair files (phased or unphased), based on the reference.
def ret_pair(wildcards):
    if 'diploid_genome' in wildcards.ref:
        # Phased setting
        return f"pairing/{wildcards.sample}.{wildcards.ref}_phased.pairs.gz"
    else:
        return f"pairing/{wildcards.sample}.{wildcards.ref}.pairs.gz"
        
# reference genomes: 
# we may obtain both a diploid_genome or an nmasked genome
REFERENCES = ['diploid_genome']

# possible phasetypes of contacts
PHASETYPES = [
    strains[0],
    strains[1],
    'unphased',
    'trans'
]

PHASEFILTER = [
    '(phase1=="0") and (phase2=="0")',
    '(phase1=="1") and (phase2=="1")',
    '(phase1==".") or (phase2==".")',
    '(phase1!=phase2) and (phase1!=".") and (phase2!=".") and (phase1!="!") and (phase2!="!")'
]

PHASEDIC = dict(map(lambda i,j : (i,j) , PHASETYPES, PHASEFILTER))


# different to bwa.snakefile
# here we skip the expensive sorting with samtools after bwa mem
rule bwa_mapping:
    input:
        fq1 = "FASTQ_fastp/{sample}_R1.fastq.gz",
        fq2 = "FASTQ_fastp/{sample}_R2.fastq.gz",
        ix = "02_bwa_index/{ref}.fa.gz.bwt"
    output:
        bam = "bwa/{sample}.{ref}.bam"
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
        bam = "bwa/{sample}.{ref}.bam",
        chr_sizes = "02_bwa_index/{ref}.chromsizes"
    output:
        pairs = "pairing/{sample}.{ref}.pairs.gz"
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
        pairs = "pairing/{sample}.diploid_genome.pairs.gz"
    output:
        pairs = "pairing/{sample}.diploid_genome_phased.pairs.gz"
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
        pairs = "pairing/{sample}.{ref}.pairs.sorted.gz"
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
        pairs = "pairing/{sample}.{ref}.pairs.sorted.gz"
    output:
        pairs = "pairing/{sample}.{ref}.pairs.dedup.gz",
        stats = "pairing/{sample}.{ref}.pairs.dedup.stats"
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
        pairs = "pairing/{sample}.diploid_genome.pairs.dedup.gz"
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
            "pairing/{sample}.{ref}.pairs.dedup.stats",
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
            -o {params.odir} \
            .
        """
