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
    shell: """
        sambamba slice -o {output.bam} {input.bam} {params.region} 2> {log}
        """

rule multiBamSummary_input:
    input:
        bams = lambda wildcards: expand("split_bam/{sample}_{part}.bam", sample=control_samples,part=wildcards.part),
        bais = lambda wildcards: expand("split_bam/{sample}_{part}.bam.bai", sample=control_samples,part=wildcards.part)
    output:
        "split_deepTools_qc/multiBamSummary/{part}.input_read_coverage.bins.npz"
    params:
        labels = " ".join(control_samples),
        blacklist = lambda wildcards: "--blackListFileName {}".format(blacklist_dict[wildcards.part]) if blacklist_dict[wildcards.part]  else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads {}".format(fragmentLength)
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
        "split_deepTools_qc/multiBamSummary/{part}.ChIP_read_coverage.bins.npz"
    params:
        labels = " ".join(chip_samples),
        blacklist = lambda wildcards: "--blackListFileName {}".format(blacklist_dict[wildcards.part]) if blacklist_dict[wildcards.part]  else "",
        read_extension = "--extendReads" if pairedEnd
                         else "--extendReads {}".format(fragmentLength)
    log:
        out = "split_deepTools_qc/logs/{part}.ChIP_multiBamSummary.out",
        err = "split_deepTools_qc/logs/{part}.ChIP_multiBamSummary.err"
    benchmark:
        "split_deepTools_qc/.benchmark/{part}.ChIP_multiBamSummary.benchmark"
    threads: 24
    conda: CONDA_SHARED_ENV
    shell: multiBamSummary_cmd
