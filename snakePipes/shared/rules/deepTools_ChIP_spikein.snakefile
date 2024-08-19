### deepTools bamCompare subtract #######################################################
part=['host','spikein']

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


if bigWigType == "subtract" or bigWigType == "both":
    rule bamCompare_subtract:
        input:
            chip_bam = "split_bam/{chip_sample}_host.bam",
            chip_bai = "split_bam/{chip_sample}_host.bam.bai",
            control_bam = "split_bam/{control_name}_host.bam",
            control_bai = "split_bam/{control_name}_host.bam.bai",
            scale_factors = "split_deepTools_qc/multiBamSummary/{part}.concatenated.scaling_factors.txt"
        output:
            "split_deepTools_ChIP/bamCompare/{chip_sample}.subtract.{control_name}.scaledBY{part}.bw"
        params:
            bwBinSize = bwBinSize,
            genome_size = genome_size,
            ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
            read_extension = "--extendReads" if pairedEnd else "--extendReads {}".format(fragmentLength),
            blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
            scaleFactors = lambda wildcards,input: " --scaleFactors {}:{} ".format(get_scaling_factor(wildcards.chip_sample,input.scale_factors),get_scaling_factor(wildcards.control_name,input.scale_factors))
        benchmark:
            "split_deepTools_ChIP/.benchmark/bamCompare.subtract.{chip_sample}.subtract.{control_name}.scaledBY{part}.benchmark"
        threads: lambda wildcards: 16 if 16<max_thread else max_thread
        conda: CONDA_SHARED_ENV
        shell: bamcompare_subtract_cmd

### deepTools bamCompare log2ratio #######################################################
if bigWigType == "log2ratio" or bigWigType == "both":
    rule bamCompare_log2:
        input:
            chip_bam = "split_bam/{chip_sample}_host.bam",
            chip_bai = "split_bam/{chip_sample}_host.bam.bai",
            control_bam = "split_bam/{control_name}_host.bam",
            control_bai = "split_bam/{control_name}_host.bam.bai",
            scale_factors = "split_deepTools_qc/multiBamSummary/{part}.concatenated.scaling_factors.txt"
        output:
            "split_deepTools_ChIP/bamCompare/{chip_sample}.log2ratio.over_{control_name}.scaledBY{part}.bw"
        params:
            bwBinSize = bwBinSize,
            ignoreForNorm = "--ignoreForNormalization {}".format(ignoreForNormalization) if ignoreForNormalization else "",
            read_extension = "--extendReads" if pairedEnd else "--extendReads {}".format(fragmentLength),
            blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
            scaleFactors = lambda wildcards,input: " --scaleFactors {}:{} ".format(get_scaling_factor(wildcards.chip_sample,input.scale_factors),get_scaling_factor(wildcards.control_name,input.scale_factors))
        benchmark:
            "split_deepTools_ChIP/.benchmark/bamCompare.log2ratio.{chip_sample}.{control_name}.scaledBY{part}.benchmark"
        threads: lambda wildcards: 16 if 16<max_thread else max_thread
        conda: CONDA_SHARED_ENV
        shell: bamcompare_log2_cmd


### deepTools plotEnrichment ###################################################

rule plotEnrichment:
    input:
        bams = expand("split_bam/{sample}_host.bam", sample = all_samples),
        bais = expand("split_bam/{sample}_host.bam.bai", sample = all_samples)
    output:
        png = "split_deepTools_ChIP/plotEnrichment/plotEnrichment.gene_features.png",
        tsv = "split_deepTools_ChIP/plotEnrichment/plotEnrichment.gene_features.tsv"
    params:
        genes_gtf = genes_gtf,
        labels = " --labels " + " ".join(all_samples),
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
        read_extension = "--extendReads" if pairedEnd else "--extendReads {}".format(fragmentLength)
    benchmark:
        "split_deepTools_ChIP/.benchmark/plotEnrichment.benchmark"
    threads: lambda wildcards: 24 if 24<max_thread else max_thread
    conda: CONDA_SHARED_ENV
    shell: plotEnrich_chip_cmd


### deepTools plotFingerprint (all files) ######################################

rule plotFingerprint:
    input:
        bams = expand("split_bam/{sample}_host.bam", sample = all_samples),
        bais = expand("split_bam/{sample}_host.bam.bai", sample = all_samples)
    output:
        metrics = "split_deepTools_ChIP/plotFingerprint/plotFingerprint.metrics.txt"
    params:
        labels = " --smartLabels ",
        blacklist = "--blackListFileName {}".format(blacklist_bed) if blacklist_bed else "",
        read_extension = "--extendReads" if pairedEnd else "--extendReads {}".format(fragmentLength),
        png = "--plotFile split_deepTools_ChIP/plotFingerprint/plotFingerprint.png" if (len(samples)<=20)
              else "",
        jsd = "--JSDsample split_bam/{}_host.bam".format(control_samples[0]) if (len(control_samples)>0)
            else ""
    benchmark:
        "split_deepTools_ChIP/.benchmark/plotFingerprint.benchmark"
    threads: lambda wildcards: 24 if 24<max_thread else max_thread
    conda: CONDA_SHARED_ENV
    shell: plotFingerprint_cmd
