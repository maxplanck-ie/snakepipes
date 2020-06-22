part=['host','spikein']
blacklist_dict={"host": blacklist_bed,"spikein": blacklist_bed_spikein}
region_dict={"host": " ".join(host_chr),"spikein": " ".join(spikein_chr)}


def get_scaling_factor(sample,input):
    sample_names=[]
    scale_factors=[]
    with open(os.path.join(outdir,input)) as f:
        for idx, line in enumerate(f):
            if idx > 0:
                sample_names.append(line.split('\t')[0])
                scale_factors.append(line.split('\t')[1])
    scale_factor = scale_factors[sample in sample_names]        

    return 1/float(scale_factor)

rule split_bamfiles_by_genome:
    input: 
        bam = "filtered_bam/{sample}.filtered.bam",
        bai = "filtered_bam/{sample}.filtered.bam.bai"
    output:
        bam = "split_bam/{sample}_{part}.bam",
        bai = "split_bam/{sample}_{part}.bam.bai"
    params:
        region = lambda wildcards: region_dict[wildcards.part]
    log: "split_bam/logs/{sample}_{part}.log"
    conda: CONDA_SAMBAMBA_ENV
    threads: 4
    shell: """
        sambamba slice -o {output.bam} {input.bam} {params.region} 2> {log};
        sambamba index -t {threads} {output.bam} 2>> {log}
        """

rule multiBamSummary_input:
    input:
        bams = lambda wildcards: expand("split_bam/{sample}_{part}.bam", sample=control_samples,part=wildcards.part),
        bais = lambda wildcards: expand("split_bam/{sample}_{part}.bam.bai", sample=control_samples,part=wildcards.part)
    output:
        npz = "split_deepTools_qc/multiBamSummary/{part}.input_read_coverage.bins.npz",
        scale_factors = "split_deepTools_qc/multiBamSummary/{part}.input.scaling_factors.txt"
    params:
        labels = " ".join(control_samples),
        blacklist = lambda wildcards: "--blackListFileName {}".format(blacklist_dict[wildcards.part]) if blacklist_dict[wildcards.part]  else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads {}".format(fragmentLength),
        scaling_factors = "--scalingFactors split_deepTools_qc/multiBamSummary/{part}.input.scaling_factors.txt",
        binSize = lambda wildcards: " --binSize 100000 " if wildcards.part=="spikein" else ""
    log:
        out = "split_deepTools_qc/logs/{part}.input_multiBamSummary.out",
        err = "split_deepTools_qc/logs/{part}.input_multiBamSummary.err"
    benchmark:
        "split_deepTools_qc/.benchmark/{part}.input_multiBamSummary.benchmark"
    threads: 24
    conda: CONDA_SHARED_ENV
    shell: multiBamSummary_cmd


rule multiBamSummary_ChIP:
    input:
        bams = lambda wildcards: expand("split_bam/{sample}_{part}.bam", sample=chip_samples,part=wildcards.part),
        bais = lambda wildcards: expand("split_bam/{sample}_{part}.bam.bai", sample=chip_samples,part=wildcards.part)
    output:
        npz = "split_deepTools_qc/multiBamSummary/{part}.ChIP_read_coverage.bins.npz",
        scale_factors = "split_deepTools_qc/multiBamSummary/{part}.ChIP.scaling_factors.txt"
    params:
        labels = " ".join(chip_samples),
        blacklist = lambda wildcards: "--blackListFileName {}".format(blacklist_dict[wildcards.part]) if blacklist_dict[wildcards.part]  else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads {}".format(fragmentLength),
        scaling_factors = "--scalingFactors split_deepTools_qc/multiBamSummary/{part}.ChIP.scaling_factors.txt",
        binSize = lambda wildcards: " --binSize 100000 " if wildcards.part=="spikein" else ""
    log:
        out = "split_deepTools_qc/logs/{part}.ChIP_multiBamSummary.out",
        err = "split_deepTools_qc/logs/{part}.ChIP_multiBamSummary.err"
    benchmark:
        "split_deepTools_qc/.benchmark/{part}.ChIP_multiBamSummary.benchmark"
    threads: 24
    conda: CONDA_SHARED_ENV
    shell: multiBamSummary_cmd


rule concatenate_scaling_factors:
    input:
        scale_factors_input = "split_deepTools_qc/multiBamSummary/{part}.input.scaling_factors.txt",
        scale_factors_chip = "split_deepTools_qc/multiBamSummary/{part}.ChIP.scaling_factors.txt"
    output: "split_deepTools_qc/multiBamSummary/{part}.concatenated.scaling_factors.txt"
    log: "split_deepTools_qc/logs/{part}.cat.scaling_factors.log"
    shell: """
        cat {input.scale_factors_input} {input.scale_factors_chip} > {output} 2> {log}
    """


rule bamCoverage_by_part:
    input:
        bam = "split_bam/{sample}_host.bam" ,
        bai = "split_bam/{sample}_host.bam.bai",
        scale_factors = "split_deepTools_qc/multiBamSummary/{part}.concatenated.scaling_factors.txt" 
    output:
        "bamCoverage/{sample}.host.seq_depth_norm.BY{part}.bw"
    params:
        bwBinSize = bwBinSize,
        genome_size = int(genome_size),
        ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads {}".format(fragmentLength),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed
                    else "",
        scaling_factors = lambda wildcards,input: "--scaleFactor {}".format(get_scaling_factor(sample,input.scale_factors)) ## subset for the one factor needed
    log:
        out = "bamCoverage/logs/bamCoverage.{sample}.BY{part}.filtered.out",
        err = "bamCoverage/logs/bamCoverage.{sample}.BY{part}.filtered.err"
    benchmark:
        "bamCoverage/.benchmark/bamCoverage.{sample}.BY{part}.filtered.benchmark"
    threads: 16  # 4GB per core
    conda: CONDA_SHARED_ENV
    shell: bamcov_cmd


rule bamPE_fragment_size:
    input:
        bams = expand("split_bam/{sample}_host.bam", sample=samples),
        bais = expand("split_bam/{sample}_host.bam.bai", sample=samples)
    output:
        "split_deepTools_qc/bamPEFragmentSize/host.fragmentSize.metric.tsv"
    params:
        plotcmd = "" if plotFormat == 'None' else
                "-o split_deepTools_qc/bamPEFragmentSize/host.fragmentSizes.{}".format(plotFormat)
    log:
        out = "split_deepTools_qc/logs/bamPEFragmentSize.out",
        err = "split_deepTools_qc/logs/bamPEFragmentSize.err"
    threads: 24
    conda: CONDA_SHARED_ENV
    shell: bamPEFragmentSize_cmd

