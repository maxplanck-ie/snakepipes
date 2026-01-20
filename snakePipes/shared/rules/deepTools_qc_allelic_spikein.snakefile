from pathlib import Path
tools_dir = Path(maindir) / "shared" / "tools"

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


#rule multiply_size_factors:
#    input:
#        allelic_sf = "deepTools_qc/multiBamSummary/allelic.scaling_factors.txt",
#        spikein_sf = "split_deepTools_qc/multiBamSummary/spikein.scaling_factors.txt"
#    output:
#        "sizeFactor_product/allelicXspikein.scaling_factors.txt"
#    params:
#        script = (tools_dir / "merge_and_multiply_and_write_sizeFactors.py"),
#        suffixes = ".genome1,.genome2",
#        a_id = "sample",
#        a_val = "scalingFactor",
#        b_id = "sample",
#        b_val = "scalingFactor",
#        float_format = ".6f"
#    script: "{params.script}"

rule bamCoverage_spikein:
    input:
        bam = "allelic_bams/{sample}.{suffix}.sorted.bam" ,
        bai = "allelic_bams/{sample}.{suffix}.sorted.bam.bai",
        scale_factors = "split_deepTools_qc/multiBamSummary/spikein.scaling_factors.txt"
    output:
        "bamCoverage/allele_specific/{sample}.{suffix}.host_scaled.BYspikein.bw"
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
        "bamCoverage/allele_specific/.benchmark/bamCoverage.{sample}.{suffix}_BYspikein.benchmark"
    threads: lambda wildcards: 16 if 16<max_thread else max_thread  # 4GB per core
    conda: CONDA_SHARED_ENV
    shell: bamcov_spikein_cmd
