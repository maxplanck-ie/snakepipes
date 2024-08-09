sample_name = os.path.splitext(os.path.basename(sampleSheet))[0]
change_direction = ["UP", "DOWN", "MIXED"]

def getInputPeaks(peakCaller, chip_samples, genrichDict):
    if peakCaller == "MACS2":
        if pipeline in 'ATACseq':
            return expand("MACS2/{chip_sample}.filtered.short.BAM_peaks.xls", chip_sample = chip_samples)
        elif pipeline == "chipseq" and useSpikeInForNorm:
            return expand("MACS2/{chip_sample}_host.BAM_peaks.xls", chip_sample = chip_samples)
        else:
            return expand("MACS2/{chip_sample}.filtered.BAM_peaks.xls", chip_sample = chip_samples)
    elif peakCaller == "HMMRATAC":
        return expand("HMMRATAC/{chip_sample}_peaks.gappedPeak", chip_sample = chip_samples)
    elif peakCaller == "SEACR":
        if pipeline == "chipseq" and useSpikeInForNorm:
            return expand("SEACR/{chip_sample}_host.stringent.bed",chip_sample=chip_samples)
        elif pipeline == "chipseq" and not useSpikeInForNorm:
            return expand("SEACR/{chip_sample}.filtered.stringent.bed",chip_sample=chip_samples)
    elif peakCaller == "Genrich":
        return expand("Genrich/{genrichGroup}.narrowPeak", genrichGroup = genrichDict.keys())
    elif externalBed:
        return externalBed


def getSizeMetrics():
    if pairedEnd:
        if not useSpikeInForNorm:
            return "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv"
        else:
            return "split_deepTools_qc/bamPEFragmentSize/host.fragmentSize.metric.tsv"
    else:
        return []

def getScaleFactors():
    if getSizeFactorsFrom=="genome":
        return "split_deepTools_qc/multiBamSummary/spikein.ChIP.scaling_factors.txt"
    elif getSizeFactorsFrom=="TSS":
        return "split_deepTools_qc/multiBamSummary_BED/spikein.ChIP.scaling_factors.txt"
    elif getSizeFactorsFrom=="input":
        return "split_deepTools_qc/multiBamSummary/spikein.input.scaling_factors.txt"
    else:
        return []

def getBamCoverage():
    if getSizeFactorsFrom=="genome":
        return expand("bamCoverage/{chip_sample}.host_scaled.BYspikein.bw", chip_sample=reordered_dict.keys())
    elif getSizeFactorsFrom=="TSS":
        return expand("bamCoverage_TSS/{chip_sample}.host_scaled.BYspikein.bw", chip_sample=reordered_dict.keys())
    elif getSizeFactorsFrom=="input":
        return expand("bamCoverage_input/{chip_sample}.host_scaled.BYspikein.bw", chip_sample=reordered_dict.keys())
    else:
        return []

def getHeatmapInput():
    if pipeline in 'ATACseq':
        return(expand("CSAW_{}_{}".format(peakCaller, sample_name) + "/CSAW.{change_dir}.cov.heatmap.png", change_dir=['UP','DOWN']))
    elif pipeline in 'chipseq':
        if chip_samples_w_ctrl:
            return(expand("CSAW_{}_{}".format(peakCaller, sample_name) + "/CSAW.{change_dir}.cov.heatmap.png", change_dir=['UP','DOWN']) + expand("CSAW_{}_{}".format(peakCaller, sample_name) + "/CSAW.{change_dir}.log2r.heatmap.png", change_dir=['UP', 'DOWN']))
        else:
          return(expand("CSAW_{}_{}".format(peakCaller, sample_name) + "/CSAW.{change_dir}.cov.heatmap.png", change_dir=['UP','DOWN']))


## CSAW for differential binding / allele-specific binding analysis
rule CSAW:
    input:
        peaks = getInputPeaks(peakCaller, chip_samples, genrichDict),
        sampleSheet = sampleSheet,
        insert_size_metrics = getSizeMetrics(),
        scale_factors = getScaleFactors() if useSpikeInForNorm else []
    output:
        "CSAW_{}_{}/CSAW.session_info.txt".format(peakCaller, sample_name),
        "CSAW_{}_{}/DiffBinding_analysis.Rdata".format(peakCaller, sample_name),
        expand("CSAW_{}_{}".format(peakCaller, sample_name) + "/Filtered.results.{change_dir}.bed", change_dir=change_direction)
    benchmark:
        "CSAW_{}_{}/.benchmark/CSAW.benchmark".format(peakCaller, sample_name)
    params:
        sampleSheet = sampleSheet,
        outdir=os.path.join(outdir, "CSAW_{}_{}".format(peakCaller, sample_name)),
        peakCaller=peakCaller,
        fdr = fdr,
        absBestLFC=absBestLFC,
        pairedEnd = pairedEnd,
        fragmentLength = fragmentLength,
        windowSize = windowSize,
        importfunc = os.path.join("shared", "rscripts", "DB_functions.R"),
        allele_info = allele_info,
        yaml_path=lambda wildcards: samples_config if pipeline in 'chipseq' else "",
        insert_size_metrics = lambda wildcards,input: os.path.join(outdir, input.insert_size_metrics) if pairedEnd else [],
        pipeline = pipeline,
        useSpikeInForNorm = useSpikeInForNorm,
        scale_factors = lambda wildcards, input: os.path.join(outdir, input.scale_factors) if input.scale_factors else "",
        externalBed = True if externalBed else False
    conda: CONDA_ATAC_ENV
    script: "../rscripts/CSAW.R"


if chip_samples_w_ctrl:
    rule calc_matrix_log2r_CSAW:
        input:
            csaw_in = "CSAW_{}_{}/CSAW.session_info.txt".format(peakCaller, sample_name),
            bigwigs = expand("split_deepTools_ChIP/bamCompare/{chip_sample}.log2ratio.over_{control_name}.scaledBYspikein.bw", zip, chip_sample=reordered_dict.keys(), control_name=reordered_dict.values()) if useSpikeInForNorm else expand("deepTools_ChIP/bamCompare/{chip_sample}.filtered.log2ratio.over_{control_name}.bw", zip, chip_sample=reordered_dict.keys(), control_name=reordered_dict.values()),
            sampleSheet = sampleSheet
        output:
            matrix = touch("CSAW_{}_{}".format(peakCaller, sample_name)+"/CSAW.{change_dir}.log2r.matrix")
        params:
            bed_in = "CSAW_{}_{}".format(peakCaller, sample_name)+"/Filtered.results.{change_dir}.bed"
        threads: 8
        conda: CONDA_SHARED_ENV
        shell: """
            if [[ -s {params.bed_in} ]]; then
                computeMatrix scale-regions -S {input.bigwigs} -R {params.bed_in} -m 1000 -b 200 -a 200 -o {output.matrix} -p {threads}
            fi
            """


    rule plot_heatmap_log2r_CSAW:
        input:
            matrix = "CSAW_{}_{}".format(peakCaller, sample_name) + "/CSAW.{change_dir}.log2r.matrix"
        output:
            image = touch("CSAW_{}_{}".format(peakCaller, sample_name) + "/CSAW.{change_dir}.log2r.heatmap.png"),
            sorted_regions = touch("CSAW_{}_{}".format(peakCaller, sample_name) + "/CSAW.{change_dir}.log2r.sortedRegions.bed")
        params:
            smpl_label=' '.join(reordered_dict.keys())
        conda: CONDA_SHARED_ENV
        shell: """
            if [[ -s {input.matrix} ]]; then
                plotHeatmap --matrixFile {input.matrix} \
                        --outFileSortedRegions {output.sorted_regions} \
                        --outFileName {output.image} \
                        --startLabel Start --endLabel End \
                        --legendLocation lower-center \
                        -x 'Scaled peak length' --labelRotation 90 \
                        --samplesLabel {params.smpl_label} --colorMap "coolwarm"
            fi
            """


rule calc_matrix_cov_CSAW:
    input:
        csaw_in = "CSAW_{}_{}/CSAW.session_info.txt".format(peakCaller, sample_name),
        bigwigs = getBamCoverage() if useSpikeInForNorm else expand("bamCoverage/{chip_sample}.filtered.seq_depth_norm.bw", chip_sample=reordered_dict.keys()),
        sampleSheet = sampleSheet
    output:
        matrix = touch("CSAW_{}_{}".format(peakCaller, sample_name) + "/CSAW.{change_dir}.cov.matrix")
    params:
        bed_in = "CSAW_{}_{}".format(peakCaller, sample_name) + "/Filtered.results.{change_dir}.bed"
    threads: 8
    conda: CONDA_SHARED_ENV
    shell: """
        touch {log.out}
        touch {log.err}
        if [[ -s {params.bed_in} ]]; then
            computeMatrix scale-regions -S {input.bigwigs} -R {params.bed_in} \
            -m 1000 -b 200 -a 200 -o {output.matrix} -p {threads}
        fi
        """


rule plot_heatmap_cov_CSAW:
    input:
        matrix = "CSAW_{}_{}".format(peakCaller, sample_name) + "/CSAW.{change_dir}.cov.matrix"
    output:
        image = touch("CSAW_{}_{}".format(peakCaller, sample_name) + "/CSAW.{change_dir}.cov.heatmap.png"),
        sorted_regions = touch("CSAW_{}_{}".format(peakCaller, sample_name) + "/CSAW.{change_dir}.cov.sortedRegions.bed")
    params:
        smpl_label=' '.join(reordered_dict.keys())
    conda: CONDA_SHARED_ENV
    shell: """
        touch {log.out}
        touch {log.err}
        if [[ -s {input.matrix} ]]; then
            plotHeatmap --matrixFile {input.matrix} \
                    --outFileSortedRegions {output.sorted_regions} \
                    --outFileName {output.image} --startLabel Start \
                    --endLabel End --legendLocation lower-center \
                    -x 'Scaled peak length' --labelRotation 90 \
                    --samplesLabel {params.smpl_label} --colorMap "coolwarm"
        fi
        """

rule CSAW_report:
    input:
        csaw_in = "CSAW_{}_{}/CSAW.session_info.txt".format(peakCaller, sample_name),
        heatmap_in = lambda wildcards: getHeatmapInput()
    output:
        outfile="CSAW_{}_{}/CSAW.Stats_report.html".format(peakCaller, sample_name)
    params:
        pipeline=pipeline,
        fdr=fdr,
        lfc=absBestLFC,
        outdir=os.path.join(outdir, "CSAW_{}_{}".format(peakCaller, sample_name)),
        sampleSheet=sampleSheet,
        useSpikeInForNorm = useSpikeInForNorm
    conda: CONDA_ATAC_ENV
    script: "../rscripts/CSAW_report.Rmd"
