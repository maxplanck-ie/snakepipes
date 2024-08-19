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



### MACS2 peak quality control #################################################

rule MACS2_peak_qc:
    input:
        bam = "split_bam/{sample}_host.bam",
        xls = "MACS2/{sample}_host.BAM_peaks.xls"
    output:
        qc = "MACS2/{sample}_host.BAM_peaks.qc.txt"
    params:
        peaks =
            lambda wildcards: "MACS2/{}_host.BAM_peaks.broadPeak".format(wildcards.sample) if is_broad(wildcards.sample)
                              else "MACS2/{}_host.BAM_peaks.narrowPeak".format(wildcards.sample),
        genome_index = genome_index
    benchmark:
        "MACS2/.benchmark/MACS2_peak_qc.{sample}_host.filtered.benchmark"
    conda: CONDA_SHARED_ENV
    shell: """
        # get the number of peaks
        peak_count=`wc -l < {params.peaks}`

        # get the number of mapped reads
        mapped_reads=`samtools view -c -F 4 {input.bam}`

        # calculate the number of alignments overlapping the peaks
        # exclude reads flagged as unmapped (unmapped reads will be reported when using -L)
        reads_in_peaks=`samtools view -c -F 4 -L {params.peaks} {input.bam}`

        # calculate Fraction of Reads In Peaks
        frip=`bc -l <<< "$reads_in_peaks/$mapped_reads"`

        # compute peak genome coverage
        peak_len=`awk '{{total+=$3-$2}}END{{print total}}' {params.peaks}`
        genome_size=`awk '{{total+=$3-$2}}END{{print total}}' {params.genome_index}`
        genomecov=`bc -l <<< "$peak_len/$genome_size"`

        # write peak-based QC metrics to output file
        printf "peak_count\tFRiP\tpeak_genome_coverage\n%d\t%5.3f\t%6.4f\n" $peak_count $frip $genomecov > {output.qc}
        """

# TODO
# add joined deepTools plotEnrichment call for all peaks and samples in one plot

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

rule SEACR_peaks_relaxed:
    input:
        chip = "filtered_bedgraph/{chip_sample}_host.fragments.bedgraph",
        control = lambda wildcards: "filtered_bedgraph/"+get_control(wildcards.chip_sample)+"_host.fragments.bedgraph" if get_control(wildcards.chip_sample)
                 else []
    output:
        "SEACR/{chip_sample}_host.relaxed.bed"
    params:
        fdr = lambda wildcards,input: fdr if not input.control else "",
        prefix = os.path.join(outdir,"SEACR/{chip_sample}_host"),
        script=os.path.join(maindir, "shared","tools/SEACR-1.3/SEACR_1.3.sh")
    conda: CONDA_SEACR_ENV
    shell: """
        bash {params.script} {input.chip} {input.control} {params.fdr} "non" "relaxed" {params.prefix}
        """


rule SEACR_peak_stringent_qc:
    input:
        bam = "split_bam/{sample}_host.bam",
        peaks = "SEACR/{sample}_host.stringent.bed"
    output:
        qc = "SEACR/{sample}_host.stringent_peaks.qc.txt"
    params:
        genome_index = genome_index
    benchmark:
        "SEACR/.benchmark/SEACR_peak_qc.{sample}_host_stringend.benchmark"
    conda: CONDA_SHARED_ENV
    shell: """
        # get the number of peaks
        peak_count=`wc -l < {input.peaks}`
        
        #get the number of mapped reads
        mapped_reads=`samtools view -c -F 4 {input.bam}`
        
        #calculate the number of alignments overlapping the peaks
        # exclude reads flagged as unmapped (unmapped reads will be reported when using -L)
        reads_in_peaks=`samtools view -c -F 4 -L {input.peaks} {input.bam}`
        
        # calculate Fraction of Reads In Peaks
        frip=`bc -l <<< "$reads_in_peaks/$mapped_reads"`
        # compute peak genome coverage
        peak_len=`awk '{{total+=$3-$2}}END{{print total}}' {input.peaks}`
        genome_size=`awk '{{total+=$3-$2}}END{{print total}}' {params.genome_index}`
        genomecov=`bc -l <<< "$peak_len/$genome_size"`
        
        # write peak-based QC metrics to output file
        printf "peak_count\tFRiP\tpeak_genome_coverage\n%d\t%5.3f\t%6.4f\n" $peak_count $frip $genomecov > {output.qc}
        """

rule SEACR_peak_relaxed_qc:
    input:
        bam = "split_bam/{sample}_host.bam",
        peaks = "SEACR/{sample}_host.relaxed.bed"
    output:
        qc = "SEACR/{sample}_host.relaxed_peaks.qc.txt"
    params:
        genome_index = genome_index
    benchmark:
        "SEACR/.benchmark/SEACR_peak_qc.{sample}_host_relaxed.benchmark"
    conda: CONDA_SHARED_ENV
    shell: """
        # get the number of peaks
        peak_count=`wc -l < {input.peaks}`
        
        #get the number of mapped reads
        mapped_reads=`samtools view -c -F 4 {input.bam}`
        
        #calculate the number of alignments overlapping the peaks
        # exclude reads flagged as unmapped (unmapped reads will be reported when using -L)
        reads_in_peaks=`samtools view -c -F 4 -L {input.peaks} {input.bam}`
        
        # calculate Fraction of Reads In Peaks
        frip=`bc -l <<< "$reads_in_peaks/$mapped_reads"`
        # compute peak genome coverage
        peak_len=`awk '{{total+=$3-$2}}END{{print total}}' {input.peaks}`
        genome_size=`awk '{{total+=$3-$2}}END{{print total}}' {params.genome_index}`
        genomecov=`bc -l <<< "$peak_len/$genome_size"`
        
        # write peak-based QC metrics to output file
        printf "peak_count\tFRiP\tpeak_genome_coverage\n%d\t%5.3f\t%6.4f\n" $peak_count $frip $genomecov > {output.qc}
        """
