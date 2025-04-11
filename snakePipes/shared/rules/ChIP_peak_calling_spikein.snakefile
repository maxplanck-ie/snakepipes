import subprocess

part=['host','spikein']

# MACS2 should be called on already filtered, e.g. duplicate-free, BAM files
# for paired-end BAM files, sambamba markdupes is fragment-based and
# therefore superior to MACS2 mate 1-based duplicate detection


### MACS2 peak calling #########################################################

if pairedEnd:
    rule writeFragmentSize:
        input: "split_deepTools_qc/bamPEFragmentSize/host.fragmentSize.metric.tsv"
        output: "MACS2/fragmentSize.metrix.tsv"


    rule MACS2:
        input:
            chip = "split_bam/{chip_sample}_host.bam",
            insert_size_metrics = "split_deepTools_qc/bamPEFragmentSize/host.fragmentSize.metric.tsv"
        output:
            peaks = "MACS2/{chip_sample}_host.BAM_peaks.xls",
            peaksPE = "MACS2/{chip_sample}_host.BAMPE_peaks.xls"
        params:
            genome_size = str(genome_size),
            broad_calling =
                lambda wildcards: "--broad" if is_broad(wildcards.chip_sample) else "",
            control_param =
                lambda wildcards: "-c split_bam/"+get_control(wildcards.chip_sample)+"_host.bam" if get_control(wildcards.chip_sample)
                else "",
            ext_size =
                lambda wildcards: " --nomodel --extsize "+get_pe_frag_length("split_bam/"+wildcards.chip_sample+"_host.bam",
                                                                            "split_deepTools_qc/bamPEFragmentSize/host.fragmentSize.metric.tsv") \
                                                                            if not cutntag else " ",
            peakCaller_options = lambda wildcards: str(peakCallerOptions or '') if not cutntag else " -p 1e-5 ",
            bampe_options = lambda wildcards: str(BAMPEPeaks or '')if not cutntag else " ",
            bam_options = lambda wildcards: str(BAMPeaks or '') if not cutntag else " "
        benchmark:
            "MACS2/.benchmark/MACS2.{chip_sample}_host.filtered.benchmark"
        conda: CONDA_CHIPSEQ_ENV
        shell: """
            macs2 callpeak -t {input.chip} {params.control_param} \
                -f BAM \
                {params.bam_options} \
                -g {params.genome_size} \
                {params.ext_size} \
                --keep-dup all \
                --outdir MACS2 \
                --name {wildcards.chip_sample}_host.BAM \
                {params.peakCaller_options} \
                {params.broad_calling}

            # also run MACS2 in paired-end mode BAMPE for comparison with single-end mode
            macs2 callpeak -t {input.chip} \
                {params.control_param} -f BAMPE \
                {params.bampe_options} \
                {params.peakCaller_options} \
                -g {params.genome_size} --keep-dup all \
                --outdir MACS2 --name {wildcards.chip_sample}_host.BAMPE \
                {params.broad_calling}
            """
else:
    rule MACS2:
        input:
            chip = "split_bam/{chip_sample}_host.bam",
            control =
                lambda wildcards: "split_bam/"+get_control(wildcards.chip_sample)+"_host.bam" if get_control(wildcards.chip_sample)
                else []
        output:
            peaks = "MACS2/{chip_sample}_host.BAM_peaks.xls",
        params:
            genome_size = int(genome_size),
            broad_calling =
                lambda wildcards: "--broad" if is_broad(wildcards.chip_sample)
                else "",
            control_param =
                lambda wildcards: "-c split_bam/"+get_control(wildcards.chip_sample)+"_host.bam" if get_control(wildcards.chip_sample)
                else "",
            frag_size=fragmentLength,
            peakCaller_options = str(peakCallerOptions or ''),
            bam_options = str(BAMPeaks or '')
        benchmark:
            "MACS2/.benchmark/MACS2.{chip_sample}_host.filtered.benchmark"
        conda: CONDA_CHIPSEQ_ENV
        shell: """
            macs2 callpeak -t {input.chip} {params.control_param} -f BAM -g {params.genome_size} \
            {params.peakCaller_options} --keep-dup all --outdir MACS2 \
            --name {wildcards.chip_sample}_host.BAM {params.bam_options} --extsize {params.frag_size} \
            {params.broad_calling}
            """


rule namesort_bams:
    input:
        bam = "split_bam/{sample}_host.bam"
    output:
        bam = temp("namesorted_bam/{sample}_host_namesorted.bam")
    params:
        tempDir = tempDir
    threads: 4
    conda: CONDA_SAMBAMBA_ENV
    shell: """
        TMPDIR={params.tempDir}
        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX)
        sambamba sort -t {threads} -o {output.bam} --tmpdir=$MYTEMP -n {input.bam}
        rm -rf $MYTEMP
         """

# Requires PE data
# Should be run once per-group!

if not isMultipleComparison:
    if pairedEnd:
        rule Genrich_peaks:
            input:
                bams=lambda wildcards: expand(os.path.join("namesorted_bam", "{sample}_host_namesorted.bam"), sample=genrichDict[wildcards.group]),
                control = lambda wildcards: ["namesorted_bam/"+get_control(x)+"_host_namesorted.bam" for x in genrichDict[wildcards.group]] if chip_samples_w_ctrl else []
            output:
                "Genrich/{group}.narrowPeak"
            params:
                bams = lambda wildcards: ",".join(expand(os.path.join("namesorted_bam", "{sample}_host_namesorted.bam"), sample=genrichDict[wildcards.group])),
                blacklist = "-E {}".format(blacklist_bed) if blacklist_bed else "",
                control_pfx=lambda wildcards,input: "-c" if input.control else "",
                control=lambda wildcards,input: ",".join(input.control) if input.control else "",
                ignoreForNorm = '-e ' + ','.join(ignoreForNormalization) if ignoreForNormalization else ""
            conda: CONDA_CHIPSEQ_ENV
            shell: """
                Genrich  -t {params.bams} {params.control_pfx} {params.control} -o {output} -r {params.blacklist} {params.ignoreForNorm} -y
                """
    else:
        rule Genrich_peaks:
            input:
                bams=lambda wildcards: expand(os.path.join("namesorted_bam", "{sample}_host_namesorted.bam"), sample=genrichDict[wildcards.group]),
                control = lambda wildcards: ["namesorted_bam/"+get_control(x)+"_host_namesorted.bam" for x in genrichDict[wildcards.group]] if chip_samples_w_ctrl else []
            output:
                "Genrich/{group}.narrowPeak"
            params:
                bams = lambda wildcards: ",".join(expand(os.path.join("namesorted_bam", "{sample}_host_namesorted.bam"), sample=genrichDict[wildcards.group])),
                blacklist = "-E {}".format(blacklist_bed) if blacklist_bed else "",
                control_pfx=lambda wildcards,input: "-c" if input.control else "",
                control=lambda wildcards,input: ",".join(input.control) if input.control else "",
                frag_size=fragmentLength,
                ignoreForNorm = "-e " + ','.join(ignoreForNormalization) if ignoreForNormalization else ""
            conda: CONDA_CHIPSEQ_ENV
            shell: """
                Genrich  -t {params.bams} {params.control_pfx} {params.control} -o {output} -r {params.blacklist} -e {params.ignoreForNorm} -w {params.frag_size}
                """
else:
    if pairedEnd:
        rule Genrich_peaks:
            input:
                bams=lambda wildcards: expand(os.path.join("namesorted_bam", "{sample}_host_namesorted.bam"), sample=genrichDict[wildcards.compGroup][wildcards.group]),
                control = lambda wildcards: ["namesorted_bam/"+get_control(x)+"_host_namesorted.bam" for x in genrichDict[wildcards.compGroup][wildcards.group]] if chip_samples_w_ctrl else []
            output:
                "Genrich/{group}.{compGroup}.narrowPeak"
            params:
                bams = lambda wildcards: ",".join(expand(os.path.join("namesorted_bam", "{sample}_host_namesorted.bam"), sample=genrichDict[wildcards.compGroup][wildcards.group])),
                blacklist = "-E {}".format(blacklist_bed) if blacklist_bed else "",
                control_pfx=lambda wildcards,input: "-c" if input.control else "",
                control=lambda wildcards,input: ",".join(input.control) if input.control else "",
                ignoreForNorm = "-e " + ','.join(ignoreForNormalization) if ignoreForNormalization else ""
            conda: CONDA_CHIPSEQ_ENV
            shell: """
                Genrich -t {params.bams} {params.control_pfx} {params.control} -o {output} -r {params.blacklist} {params.ignoreForNorm} -y
                """
    else:
        rule Genrich_peaks:
            input:
                bams=lambda wildcards: expand(os.path.join("namesorted_bam", "{sample}_host_namesorted.bam"), sample=genrichDict[wildcards.compGroup][wildcards.group]),
                control = lambda wildcards: ["namesorted_bam/"+get_control(x)+"_host_namesorted.bam" for x in genrichDict[wildcards.compGroup][wildcards.group] ] if chip_samples_w_ctrl else []
            output:
                "Genrich/{group}.{compGroup}.narrowPeak"
            params:
                bams = lambda wildcards: ",".join(expand(os.path.join("namesorted_bam", "{sample}_host_namesorted.bam"), sample=genrichDict[wildcards.compGroup][wildcards.group])),
                blacklist = "-E {}".format(blacklist_bed) if blacklist_bed else "",
                control_pfx=lambda wildcards,input: "-c" if input.control else "",
                control=lambda wildcards,input: ",".join(input.control) if input.control else "",
                frag_size=fragmentLength,
                ignoreForNorm = "-e " + ','.join(ignoreForNormalization) if ignoreForNormalization else ""
            conda: CONDA_CHIPSEQ_ENV
            shell: """
                Genrich -t {params.bams} {params.control_pfx} {params.control} -o {output} -r {params.blacklist} -e {params.ignoreForNorm} -w {params.frag_size}
                """


rule prep_bedgraph:
    input: "bamCoverage/{sample}.host_scaled.BYhost.bw"
    output: temp("filtered_bedgraph/{sample}_host.fragments.bedgraph")
    conda: CONDA_SEACR_ENV
    shell: """
        bigWigToBedGraph {input} {output}
        """

rule SEACR_peaks_stringent:
    input:
        chip = "filtered_bedgraph/{chip_sample}_host.fragments.bedgraph",
        control = lambda wildcards: "filtered_bedgraph/"+get_control(wildcards.chip_sample)+"_host.fragments.bedgraph" if get_control(wildcards.chip_sample)
                 else []
    output:
        "SEACR/{chip_sample}_host.stringent.bed"
    params:
        fdr = lambda wildcards,input: fdr if not input.control else "",
        prefix = os.path.join(outdir,"SEACR/{chip_sample}_host"),
        script=os.path.join(maindir, "shared","tools/SEACR-1.3/SEACR_1.3.sh")
    conda: CONDA_SEACR_ENV
    shell: """
        bash {params.script} {input.chip} {input.control} {params.fdr} "non" "stringent" {params.prefix}
        """


def collectPeaks(caller):
    if caller == "SEACR":
        return expand("SEACR/{chip_sample}_host.stringent.bed", chip_sample=chip_samples)
    elif caller == "MACS2":
        return expand("MACS2/{chip_sample}_host.BAM_peaks.xls",chip_sample=chip_samples)
    elif caller == "Genrich":
        return expand("Genrich/{group}.narrowPeak",group=genrichDict.keys())


rule chipqc:
    input:
        bams = expand("split_bam/{chip_sample}_host.bam",chip_sample=chip_samples),
        peaks = collectPeaks(caller=peakCaller),
        sampleSheet = sampleSheet if sampleSheet else [],
        chipdict = os.path.join(outdir,"chip_samples.yaml")
    output:
        "{}_chipqc/sessionInfo.txt".format(peakCaller)
    params:
        genome = genome,
        outdir = "{}_chipqc".format(peakCaller),
        blacklist = blacklist_bed,
        bams = lambda wildcards,input: [os.path.join(outdir,x) for x in input.bams],
        peaks = lambda wildcards,input: [os.path.join(outdir,x) for x in input.peaks],
        narrow_samples = narrow_samples,
        broad_samples = broad_samples,
        useSpikeinForNorm = useSpikeInForNorm
    threads: 8
    benchmark:
        "{}_chipqc/.benchmark/chipqc.benchmark".format(peakCaller)
    conda: CONDA_CHIPQC_ENV
    script: "../rscripts/chipqc.R"
