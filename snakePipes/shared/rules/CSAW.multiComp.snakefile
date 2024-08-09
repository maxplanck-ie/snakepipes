#sample_name = os.path.splitext(os.path.basename(sampleSheet))[0]
change_direction = ["UP", "DOWN", "MIXED"]
#compGroup = cf.returnComparisonGroups(sampleSheet)

def get_outdir(peak_caller,sampleSheet):
    sample_name = os.path.splitext(os.path.basename(str(sampleSheet)))[0]

    return("CSAW_{}_{}".format(peak_caller, sample_name))

def getInputPeaks(peakCaller, chip_samples, genrichDict,comp_group):
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
        return expand("Genrich/{genrichGroup}.{{compGroup}}.narrowPeak", genrichGroup = genrichDict[comp_group].keys())
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

def getBamCoverage(comp_group):
    if getSizeFactorsFrom=="genome":
        return expand("bamCoverage/{chip_sample}.host_scaled.BYspikein.bw", chip_sample=reordered_dict[comp_group].keys())
    elif getSizeFactorsFrom=="TSS":
        return expand("bamCoverage_TSS/{chip_sample}.host_scaled.BYspikein.bw", chip_sample=reordered_dict[comp_group].keys())
    elif getSizeFactorsFrom=="input":
        return expand("bamCoverage_input/{chip_sample}.host_scaled.BYspikein.bw", chip_sample=reordered_dict[comp_group].keys())
    else:
        return []

def getHeatmapInput():
    if pipeline in 'ATACseq':
        return(expand("CSAW_{}_{}".format(peakCaller, sample_name + ".{{compGroup}}") + "/CSAW.{change_dir}.cov.heatmap.png", change_dir=['UP','DOWN']))
    elif pipeline in 'chipseq':
        if chip_samples_w_ctrl:
            return(expand("CSAW_{}_{}".format(peakCaller, sample_name + ".{{compGroup}}") + "/CSAW.{change_dir}.cov.heatmap.png", change_dir=['UP','DOWN']) + expand("CSAW_{}_{}".format(peakCaller, sample_name + ".{{compGroup}}") + "/CSAW.{change_dir}.log2r.heatmap.png", change_dir=['UP', 'DOWN']))
        else:
          return(expand("CSAW_{}_{}".format(peakCaller, sample_name + ".{{compGroup}}") + "/CSAW.{change_dir}.cov.heatmap.png", change_dir=['UP','DOWN']))

checkpoint split_sampleSheet:
    input:
        sampleSheet = sampleSheet
    output:
        splitSheets = os.path.join("splitSampleSheets",os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv")
    params:
        splitSheetPfx = os.path.join("splitSampleSheets",os.path.splitext(os.path.basename(str(sampleSheet)))[0])
    run:
        if isMultipleComparison:
            cf.splitSampleSheet(input.sampleSheet,params.splitSheetPfx)

## CSAW for differential binding / allele-specific binding analysis
rule CSAW:
    input:
        peaks = lambda wildcards: getInputPeaks(peakCaller, chip_samples, genrichDict, comp_group=wildcards.compGroup),
        sampleSheet = lambda wildcards: checkpoints.split_sampleSheet.get(compGroup=wildcards.compGroup).output,
        insert_size_metrics = getSizeMetrics(),
        scale_factors = getScaleFactors() if useSpikeInForNorm else []
    output:
        "{}/CSAW.session_info.txt".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv")),
        "{}/DiffBinding_analysis.Rdata".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv")),
        expand("{}".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{{compGroup}}.tsv")) + "/Filtered.results.{change_dir}.bed", change_dir=change_direction)
    benchmark:
        "{}/.benchmark/CSAW.benchmark".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv"))
    params:
        sampleSheet = lambda wildcards,input: os.path.join(outdir, str(input.sampleSheet)),
        outdir=lambda wildcards,input: get_outdir(peakCaller,input.sampleSheet),
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


rule calc_matrix_log2r_CSAW:
    input:
        csaw_in = "{}/CSAW.session_info.txt".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv")),
        bigwigs = lambda wildcards: expand("split_deepTools_ChIP/bamCompare/{chip_sample}.log2ratio.over_{control_name}.scaledBYspikein.bw", zip, chip_sample=reordered_dict[wildcards.compGroup].keys(), control_name=reordered_dict[wildcards.compGroup].values()) if useSpikeInForNorm else expand("deepTools_ChIP/bamCompare/{chip_sample}.filtered.log2ratio.over_{control_name}.bw", zip, chip_sample=reordered_dict[wildcards.compGroup].keys(), control_name=reordered_dict[wildcards.compGroup].values()),
        sampleSheet = sampleSheet
    output:
        matrix = touch("{}".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv"))+"/CSAW.{change_dir}.log2r.matrix")
    params:
        bed_in = "{}".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv"))+"/Filtered.results.{change_dir}.bed"
    threads: 8
    conda: CONDA_SHARED_ENV
    shell: """
        touch {log.out}
        touch {log.err}
        if [[ -s {params.bed_in} ]]; then
            computeMatrix scale-regions -S {input.bigwigs} -R {params.bed_in} -m 1000 -b 200 -a 200 -o {output.matrix} -p {threads}
        fi
        """


rule plot_heatmap_log2r_CSAW:
    input:
        matrix = "{}".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv")) + "/CSAW.{change_dir}.log2r.matrix"
    output:
        image = touch("{}".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv")) + "/CSAW.{change_dir}.log2r.heatmap.png"),
        sorted_regions = touch("{}".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv")) + "/CSAW.{change_dir}.log2r.sortedRegions.bed")
    params:
        smpl_label = lambda wildcards: ' '.join(reordered_dict[wildcards.compGroup].keys())
    conda: CONDA_SHARED_ENV
    shell: """
        touch {log.out}
        touch {log.err}
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
        csaw_in = "{}/CSAW.session_info.txt".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv")),
        bigwigs = lambda wildcards: getBamCoverage(wildcards.compGroup) if useSpikeInForNorm else expand("bamCoverage/{chip_sample}.filtered.seq_depth_norm.bw", chip_sample=reordered_dict[wildcards.compGroup].keys()),
        sampleSheet = sampleSheet
    output:
        matrix = touch("{}".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv")) + "/CSAW.{change_dir}.cov.matrix")
    params:
        bed_in = "{}".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv")) + "/Filtered.results.{change_dir}.bed"
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
        matrix = "{}".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv")) + "/CSAW.{change_dir}.cov.matrix"
    output:
        image = touch("{}".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv")) + "/CSAW.{change_dir}.cov.heatmap.png"),
        sorted_regions = touch("{}".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv")) + "/CSAW.{change_dir}.cov.sortedRegions.bed")
    params:
        smpl_label = lambda wildcards: ' '.join(reordered_dict[wildcards.compGroup].keys())
    conda: CONDA_SHARED_ENV
    shell: """
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
        csaw_in = "{}/CSAW.session_info.txt".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv")),
        heatmap_in = lambda wildcards: getHeatmapInput()
    output:
        outfile="{}/CSAW.Stats_report.html".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv"))
    params:
        pipeline=pipeline,
        fdr=fdr,
        lfc=absBestLFC,
        outdir=os.path.join(outdir, "{}".format(get_outdir(peakCaller,os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv"))),
        sampleSheet=sampleSheet,
        useSpikeInForNorm = useSpikeInForNorm
    conda: CONDA_ATAC_ENV
    script: "../rscripts/CSAW_report.Rmd"
