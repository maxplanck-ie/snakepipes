part=['host','spikein']
blacklist_dict={"host": blacklist_bed,"spikein": spikein_blacklist_bed }
region_dict={"host": " ".join(host_chr.keys()),"spikein": " ".join(spikein_chr.keys())}


def get_scaling_factor(sample,input):
    sample_names=[]
    scale_factors=[]
    if os.path.isfile(os.path.join(outdir,input)):
        with open(os.path.join(outdir,input)) as f:
            for idx, line in enumerate(f):
                if idx > 0:
                    sample_names.append(line.split('\t')[0])
                    scale_factors.append((line.split('\t')[1]).rstrip("\n"))
        sf_dict = dict(zip(sample_names, scale_factors))
        scale_factor = sf_dict[sample]

        return float(scale_factor)
    else:
        return float(1)

rule split_bamfiles_by_genome:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        bai = "filtered_bam/{sample}.filtered.bam.bai"
    output:
        bam = "split_bam/{sample}_{part}.bam",
        bai = "split_bam/{sample}_{part}.bam.bai"
    params:
        region = lambda wildcards: region_dict[wildcards.part]
    conda: CONDA_SAMBAMBA_ENV
    threads: 4
    shell: """
        sambamba slice -o {output.bam} {input.bam} {params.region};
        sambamba index -t {threads} {output.bam}
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
        binSize = lambda wildcards: " --binSize "+str(spikein_bin_size) if wildcards.part=="spikein" else "",
        spikein_region = lambda wildcards: " --region "+spikein_region if ((wildcards.part=="spikein") and (spikein_region != "")) else ""
    benchmark:
        "split_deepTools_qc/.benchmark/{part}.input_multiBamSummary.benchmark"
    threads: lambda wildcards: 24 if 24<max_thread else max_thread
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
        binSize = lambda wildcards: " --binSize "+str(spikein_bin_size) if wildcards.part=="spikein" else "",
        spikein_region = lambda wildcards: " --region "+spikein_region if ((wildcards.part=="spikein") and (spikein_region != "")) else ""
    benchmark:
        "split_deepTools_qc/.benchmark/{part}.ChIP_multiBamSummary.benchmark"
    threads: lambda wildcards: 24 if 24<max_thread else max_thread
    conda: CONDA_SHARED_ENV
    shell: multiBamSummary_cmd


rule multiBamSummary_TSS:
    input:
        bams = lambda wildcards: expand("split_bam/{sample}_spikein.bam", sample=chip_samples),
        bais = lambda wildcards: expand("split_bam/{sample}_spikein.bam.bai", sample=chip_samples),
        bed = "Annotation_spikein/TSS.filtered.bed"
    output:
        npz = "split_deepTools_qc/multiBamSummary_BED/spikein.ChIP_read_coverage.bins.npz",
        scale_factors = "split_deepTools_qc/multiBamSummary_BED/spikein.ChIP.scaling_factors.txt"
    params:
        labels = " ".join(chip_samples),
        blacklist = lambda wildcards: "--blackListFileName {}".format(blacklist_dict['spikein']) if blacklist_dict['spikein']  else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads {}".format(fragmentLength),
        scaling_factors = "--scalingFactors split_deepTools_qc/multiBamSummary_BED/spikein.ChIP.scaling_factors.txt",
        binSize = " --binSize 100000 "
    benchmark:
        "split_deepTools_qc/.benchmark/spikein.ChIP_multiBamSummary.BED.benchmark"
    threads: lambda wildcards: 24 if 24<max_thread else max_thread
    conda: CONDA_SHARED_ENV
    shell: multiBamSummary_spikein_cmd


rule concatenate_scaling_factors:
    input:
        scale_factors_input = "split_deepTools_qc/multiBamSummary/{part}.input.scaling_factors.txt",
        scale_factors_chip = "split_deepTools_qc/multiBamSummary/{part}.ChIP.scaling_factors.txt"
    output: "split_deepTools_qc/multiBamSummary/{part}.concatenated.scaling_factors.txt"
    shell: """
        cat {input.scale_factors_input} {input.scale_factors_chip} > {output}
    """


rule bamCoverage_by_part:
    input:
        bam = "split_bam/{sample}_host.bam" ,
        bai = "split_bam/{sample}_host.bam.bai",
        scale_factors = "split_deepTools_qc/multiBamSummary/{part}.ChIP.scaling_factors.txt"
    output:
        "bamCoverage/{sample}.host_scaled.BY{part}.bw"
    params:
        bwBinSize = bwBinSize,
        genome_size = int(genome_size),
        ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads {}".format(fragmentLength),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed
                    else "",
        scaling_factors = lambda wildcards,input: "--scaleFactor {}".format(get_scaling_factor(wildcards.sample,input.scale_factors)) ## subset for the one factor needed
    benchmark:
        "bamCoverage/.benchmark/bamCoverage.{sample}.BY{part}.filtered.benchmark"
    threads: lambda wildcards: 16 if 16<max_thread else max_thread  # 4GB per core
    conda: CONDA_SHARED_ENV
    shell: bamcov_spikein_cmd


rule bamCoverage_by_TSS:
    input:
        bam = "split_bam/{sample}_host.bam" ,
        bai = "split_bam/{sample}_host.bam.bai",
        scale_factors = "split_deepTools_qc/multiBamSummary_BED/spikein.ChIP.scaling_factors.txt"
    output:
        "bamCoverage_TSS/{sample}.host_scaled.BYspikein.bw"
    params:
        bwBinSize = bwBinSize,
        genome_size = int(genome_size),
        ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads {}".format(fragmentLength),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed
                    else "",
        scaling_factors = lambda wildcards,input: "--scaleFactor {}".format(get_scaling_factor(wildcards.sample,input.scale_factors)) ## subset for the one factor needed
    benchmark:
        "bamCoverage_TSS/.benchmark/bamCoverage.{sample}.BYspikein.filtered.benchmark"
    threads: lambda wildcards: 16 if 16<max_thread else max_thread  # 4GB per core
    conda: CONDA_SHARED_ENV
    shell: bamcov_spikein_cmd


rule bamCoverage_by_input:
    input:
        bam = "split_bam/{sample}_host.bam" ,
        bai = "split_bam/{sample}_host.bam.bai",
        scale_factors = "split_deepTools_qc/multiBamSummary/spikein.input.scaling_factors.txt"
    output:
        "bamCoverage_input/{sample}.host_scaled.BYspikein.bw"
    params:
        bwBinSize = bwBinSize,
        genome_size = int(genome_size),
        ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads {}".format(fragmentLength),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed
                    else "",
        scaling_factors = lambda wildcards,input: "--scaleFactor {}".format(get_scaling_factor(get_control(wildcards.sample),input.scale_factors)) ## subset for the one factor needed
    benchmark:
        "bamCoverage_input/.benchmark/bamCoverage.{sample}.BYspikein.filtered.benchmark"
    threads: lambda wildcards: 16 if 16<max_thread else max_thread  # 4GB per core
    conda: CONDA_SHARED_ENV
    shell: bamcov_spikein_cmd


rule bamPE_fragment_size:
    input:
        bams = lambda wildcards: expand("split_bam/{sample}_{part}.bam", sample=samples,part=wildcards.part),
        bais = lambda wildcards: expand("split_bam/{sample}_{part}.bam.bai", sample=samples,part=wildcards.part)
    output:
        "split_deepTools_qc/bamPEFragmentSize/{part}.fragmentSize.metric.tsv"
    params:
        plotcmd = lambda wildcards: "" if plotFormat == 'None' else
                "-o split_deepTools_qc/bamPEFragmentSize/" + wildcards.part + ".fragmentSizes.{}".format(plotFormat)
    threads: lambda wildcards: 24 if 24<max_thread else max_thread
    conda: CONDA_SHARED_ENV
    shell: bamPEFragmentSize_cmd
