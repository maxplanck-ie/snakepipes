part=['host','spikein']
blacklist_dict={"host": blacklist_bed,"spikein": blacklist_bed_spikein}
region_dict={"host": " ".join(host_chr),"spikein": " ".join(spikein_chr)}

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
        scaling_factors = "--scalingFactors split_deepTools_qc/multiBamSummary/{part}.input.scaling_factors.txt"
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
        scaling_factors = "--scalingFactors split_deepTools_qc/multiBamSummary/{part}.ChIP.scaling_factors.txt"
    log:
        out = "split_deepTools_qc/logs/{part}.ChIP_multiBamSummary.out",
        err = "split_deepTools_qc/logs/{part}.ChIP_multiBamSummary.err"
    benchmark:
        "split_deepTools_qc/.benchmark/{part}.ChIP_multiBamSummary.benchmark"
    threads: 24
    conda: CONDA_SHARED_ENV
    shell: multiBamSummary_cmd

rule bamCoverage_input_by_host:
    input:
        bams = lambda wildcards: expand("split_bam/{sample}_host.bam", sample=control_samples),
        bais = lambda wildcards: expand("split_bam/{sample}_host.bam.bai", sample=control_samples),
        scale_factors = "split_deepTools_qc/multiBamSummary/host.input.scaling_factors.txt"
    output:
        "bamCoverage_NormedByHost/{sample}.host.seq_depth_norm.bw"
    params:
        bwBinSize = bwBinSize,
        genome_size = int(genome_size),
        ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads {}".format(fragmentLength),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed
                    else "",
        scaling_factors = "--scaleFactor split_deepTools_qc/multiBamSummary/host.input.scaling_factors.txt"
    log:
        out = "bamCoverage_NormedByHost/logs/bamCoverage.{sample}.filtered.out",
        err = "bamCoverage_NormedByHost/logs/bamCoverage.{sample}.filtered.err"
    benchmark:
        "bamCoverage_NormedByHost/.benchmark/bamCoverage.{sample}.filtered.benchmark"
    threads: 16  # 4GB per core
    conda: CONDA_SHARED_ENV
    shell: bamcov_cmd

rule bamCoverage_input_by_spikein:
    input:
        bams = lambda wildcards: expand("split_bam/{sample}_host.bam", sample=control_samples),
        bais = lambda wildcards: expand("split_bam/{sample}_host.bam.bai", sample=control_samples),
        scale_factors = "split_deepTools_qc/multiBamSummary/spikein.input.scaling_factors.txt"
    output:
        "bamCoverage_NormedBySpikeIn/{sample}.spikein.seq_depth_norm.bw"
    params:
        bwBinSize = bwBinSize,
        genome_size = int(genome_size),
        ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads {}".format(fragmentLength),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed
                    else "",
        scaling_factors = "--scaleFactor split_deepTools_qc/multiBamSummary/spikein.input.scaling_factors.txt"
    log:
        out = "bamCoverage_NormedBySpikeIn/logs/bamCoverage.{sample}.filtered.out",
        err = "bamCoverage_NormedBySpikeIn/logs/bamCoverage.{sample}.filtered.err"
    benchmark:
        "bamCoverage_NormedBySpikein/.benchmark/bamCoverage.{sample}.filtered.benchmark"
    threads: 16  # 4GB per core
    conda: CONDA_SHARED_ENV
    shell: bamcov_cmd

rule bamCoverage_ChIP_by_host:
    input:
        bams = lambda wildcards: expand("split_bam/{sample}_host.bam", sample=chip_samples),
        bais = lambda wildcards: expand("split_bam/{sample}_host.bam.bai", sample=chip_samples),
        scale_factors = "split_deepTools_qc/multiBamSummary/host.ChIP.scaling_factors.txt"
    output:
        "bamCoverage_NormedByHost/{sample}.host.seq_depth_norm.bw"
    params:
        bwBinSize = bwBinSize,
        genome_size = int(genome_size),
        ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads {}".format(fragmentLength),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed
                    else "",
        scaling_factors = "--scaleFactor split_deepTools_qc/multiBamSummary/host.input.scaling_factors.txt"
    log:
        out = "bamCoverage_NormedByHost/logs/bamCoverage.{sample}.filtered.out",
        err = "bamCoverage_NormedByHost/logs/bamCoverage.{sample}.filtered.err"
    benchmark:
        "bamCoverage_NormedByHost/.benchmark/bamCoverage.{sample}.filtered.benchmark"
    threads: 16  # 4GB per core
    conda: CONDA_SHARED_ENV
    shell: bamcov_cmd

rule bamCoverage_ChIP_by_spikein:
    input:
        bams = lambda wildcards: expand("split_bam/{sample}_host.bam", sample=chip_samples),
        bais = lambda wildcards: expand("split_bam/{sample}_host.bam.bai", sample=chip_samples),
        scale_factors = "split_deepTools_qc/multiBamSummary/spikein.ChIP.scaling_factors.txt"
    output:
        "bamCoverage_NormedBySpikeIn/{sample}.spikein.seq_depth_norm.bw"
    params:
        bwBinSize = bwBinSize,
        genome_size = int(genome_size),
        ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads {}".format(fragmentLength),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed
                    else "",
        scaling_factors = "--scaleFactor split_deepTools_qc/multiBamSummary/spikein.ChIP.scaling_factors.txt"
    log:
        out = "bamCoverage_NormedBySpikeIn/logs/bamCoverage.{sample}.filtered.out",
        err = "bamCoverage_NormedBySpikeIn/logs/bamCoverage.{sample}.filtered.err"
    benchmark:
        "bamCoverage_NormedBySpikein/.benchmark/bamCoverage.{sample}.filtered.benchmark"
    threads: 16  # 4GB per core
    conda: CONDA_SHARED_ENV
    shell: bamcov_cmd
