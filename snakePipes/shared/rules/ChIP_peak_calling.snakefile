import subprocess

# MACS2 should be called on already filtered, e.g. duplicate-free, BAM files
# for paired-end BAM files, sambamba markdupes is fragment-based and
# therefore superior to MACS2 mate 1-based duplicate detection


### MACS2 peak calling #########################################################

if pairedEnd:
    rule writeFragmentSize:
        input: "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv"
        output: "MACS2/fragmentSize.metrix.tsv"


    rule MACS2:
        input:
            chip = "filtered_bam/{chip_sample}.filtered.bam",
            frag_size = "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv"
        output:
            peaks = "MACS2/{chip_sample}.filtered.BAM_peaks.xls",
            peaksPE = "MACS2/{chip_sample}.filtered.BAMPE_peaks.xls"
        params:
            broad_calling =
                lambda wildcards: "--broad " if is_broad(wildcards.chip_sample) else "",
            control_param =
                lambda wildcards: " -c filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam" if get_control(wildcards.chip_sample)
                else "",
            genome_size = str(genome_size),
            ext_size =
                lambda wildcards: " --nomodel --extsize "+get_pe_frag_length("filtered_bam/"+wildcards.chip_sample+".filtered.bam",
                                                                            "deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv") \
                                                                            if not cutntag else " ",
            peakCaller_options = lambda wildcards: str(peakCallerOptions or '') if not cutntag else " -p 1e-5 ",
            bampe_options = lambda wildcards: str(BAMPEPeaks or '')if not cutntag else " ",
            bam_options = lambda wildcards: str(BAMPeaks or '') if not cutntag else " "
        log:
            out = "MACS2/logs/MACS2.{chip_sample}.filtered.out",
            err = "MACS2/logs/MACS2.{chip_sample}.filtered.err"
        benchmark:
            "MACS2/.benchmark/MACS2.{chip_sample}.filtered.benchmark"
        conda: CONDA_CHIPSEQ_ENV
        shell: """
            macs2 callpeak -t {input.chip} {params.control_param} \
                -f BAM \
                {params.bam_options} \
                -g {params.genome_size} \
                {params.ext_size} \
                --keep-dup all \
                --outdir MACS2 \
                --name {wildcards.chip_sample}.filtered.BAM \
                {params.peakCaller_options} \
                {params.broad_calling} > {log.out} 2> {log.err}

            # also run MACS2 in paired-end mode BAMPE for comparison with single-end mode
            macs2 callpeak -t {input.chip} \
                {params.control_param} -f BAMPE \
                {params.bampe_options} \
                {params.peakCaller_options} \
                -g {params.genome_size} --keep-dup all \
                --outdir MACS2 --name {wildcards.chip_sample}.filtered.BAMPE \
                {params.broad_calling} > {log.out}.BAMPE 2> {log.err}.BAMPE
            """
else:
    rule MACS2:
        input:
            chip = "filtered_bam/{chip_sample}.filtered.bam",
            control =
                lambda wildcards: "filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam" if get_control(wildcards.chip_sample)
                else []
        output:
            peaks = "MACS2/{chip_sample}.filtered.BAM_peaks.xls",
        params:
            genome_size = str(genome_size),
            broad_calling =
                lambda wildcards: "--broad" if is_broad(wildcards.chip_sample)
                else "",
            control_param =
                lambda wildcards: " -c filtered_bam/"+get_control(wildcards.chip_sample)+".filtered.bam" if get_control(wildcards.chip_sample)
                else "",
            frag_size=fragmentLength,
            peakCaller_options = str(peakCallerOptions or ''),
            bam_options = str(BAMPeaks or '')
        log:
            out = "MACS2/logs/MACS2.{chip_sample}.filtered.out",
            err = "MACS2/logs/MACS2.{chip_sample}.filtered.err"
        benchmark:
            "MACS2/.benchmark/MACS2.{chip_sample}.filtered.benchmark"
        conda: CONDA_CHIPSEQ_ENV
        shell: """
            macs2 callpeak -t {input.chip} {params.control_param} -f BAM -g {params.genome_size} \
            {params.peakCaller_options} --keep-dup all --outdir MACS2 \
            --name {wildcards.chip_sample}.filtered.BAM {params.bam_options} --extsize {params.frag_size} \
            {params.broad_calling} > {log.out} 2> {log.err}
            """



### MACS2 peak quality control #################################################

rule MACS2_peak_qc:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        xls = "MACS2/{sample}.filtered.BAM_peaks.xls"
    output:
        qc = "MACS2/{sample}.filtered.BAM_peaks.qc.txt"
    params:
        peaks =
            lambda wildcards: "MACS2/{}.filtered.BAM_peaks.broadPeak".format(wildcards.sample) if is_broad(wildcards.sample)
                              else "MACS2/{}.filtered.BAM_peaks.narrowPeak".format(wildcards.sample),
        genome_index = genome_index
    benchmark:
        "MACS2/.benchmark/MACS2_peak_qc.{sample}.filtered.benchmark"
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
        bam = "filtered_bam/{sample}.filtered.bam"
    output:
        bam = temp("filtered_bam/{sample}.namesorted.bam")
    log:
        "filtered_bam/logs/{sample}.namesort.err"
    params:
        tempDir = tempDir
    threads: 4
    conda: CONDA_SAMBAMBA_ENV
    shell: """
        TMPDIR={params.tempDir}
        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX)
        sambamba sort -t {threads} -o {output.bam} --tmpdir=$MYTEMP -n {input.bam} 2> {log}
        rm -rf $MYTEMP
         """

# Requires PE data
# Should be run once per-group!
if not isMultipleComparison:
    if pairedEnd:
        rule Genrich_peaks:
            input:
                bams=lambda wildcards: expand(os.path.join("filtered_bam", "{sample}.namesorted.bam"), sample=genrichDict[wildcards.group]),
                control = lambda wildcards: ["filtered_bam/"+get_control(x)+".namesorted.bam" for x in genrichDict[wildcards.group]] if chip_samples_w_ctrl else []
            output:
                "Genrich/{group}.narrowPeak"
            log: "Genrich/logs/{group}.log"
            params:
                bams = lambda wildcards: ",".join(expand(os.path.join("filtered_bam", "{sample}.namesorted.bam"), sample=genrichDict[wildcards.group])),
                blacklist = "-E {}".format(blacklist_bed) if blacklist_bed else "",
                control_pfx=lambda wildcards,input: "-c" if input.control else "",
                control=lambda wildcards,input: ",".join(input.control) if input.control else "",
                ignoreForNorm = "-e " + ','.join(ignoreForNormalization) if ignoreForNormalization else ""
            conda: CONDA_CHIPSEQ_ENV
            shell: """
                Genrich -t {params.bams} {params.control_pfx} {params.control} -o {output} -r {params.blacklist} {params.ignoreForNorm} -y 2> {log}
                """
    else:
        rule Genrich_peaks:
            input:
                bams=lambda wildcards: expand(os.path.join("filtered_bam", "{sample}.namesorted.bam"), sample=genrichDict[wildcards.group]),
                control = lambda wildcards: ["filtered_bam/"+get_control(x)+".namesorted.bam" for x in genrichDict[wildcards.group] ] if chip_samples_w_ctrl else []
            output:
                "Genrich/{group}.narrowPeak"
            log: "Genrich/logs/{group}.log"
            params:
                bams = lambda wildcards: ",".join(expand(os.path.join("filtered_bam", "{sample}.namesorted.bam"), sample=genrichDict[wildcards.group])),
                blacklist = "-E {}".format(blacklist_bed) if blacklist_bed else "",
                control_pfx=lambda wildcards,input: "-c" if input.control else "",
                control=lambda wildcards,input: ",".join(input.control) if input.control else "",
                frag_size=fragmentLength,
                ignoreForNorm = "-e " + ','.join(ignoreForNormalization) if ignoreForNormalization else ""
            conda: CONDA_CHIPSEQ_ENV
            shell: """
                Genrich -t {params.bams} {params.control_pfx} {params.control} -o {output} -r {params.blacklist} {params.ignoreForNorm} -w {params.frag_size} 2> {log}
                """
else:
    if pairedEnd:
        rule Genrich_peaks:
            input:
                bams=lambda wildcards: expand(os.path.join("filtered_bam", "{sample}.namesorted.bam"), sample=genrichDict[wildcards.compGroup][wildcards.group]),
                control = lambda wildcards: ["filtered_bam/"+get_control(x)+".namesorted.bam" for x in genrichDict[wildcards.compGroup][wildcards.group]] if chip_samples_w_ctrl else []
            output:
                "Genrich/{group}.{compGroup}.narrowPeak"
            log: "Genrich/logs/{group}.{compGroup}.log"
            params:
                bams = lambda wildcards: ",".join(expand(os.path.join("filtered_bam", "{sample}.namesorted.bam"), sample=genrichDict[wildcards.compGroup][wildcards.group])),
                blacklist = "-E {}".format(blacklist_bed) if blacklist_bed else "",
                control_pfx=lambda wildcards,input: "-c" if input.control else "",
                control=lambda wildcards,input: ",".join(input.control) if input.control else "",
                ignoreForNorm = "-e " + ','.join(ignoreForNormalization) if ignoreForNormalization else ""
            conda: CONDA_CHIPSEQ_ENV
            shell: """
                Genrich -t {params.bams} {params.control_pfx} {params.control} -o {output} -r {params.blacklist} {params.ignoreForNorm} -y 2> {log}
                """
    else:
        rule Genrich_peaks:
            input:
                bams=lambda wildcards: expand(os.path.join("filtered_bam", "{sample}.namesorted.bam"), sample=genrichDict[wildcards.compGroup][wildcards.group]),
                control = lambda wildcards: ["filtered_bam/"+get_control(x)+".namesorted.bam" for x in genrichDict[wildcards.compGroup][wildcards.group] ] if chip_samples_w_ctrl else []
            output:
                "Genrich/{group}.{compGroup}.narrowPeak"
            log: "Genrich/logs/{group}.{compGroup}.log"
            params:
                bams = lambda wildcards: ",".join(expand(os.path.join("filtered_bam", "{sample}.namesorted.bam"), sample=genrichDict[wildcards.compGroup][wildcards.group])),
                blacklist = "-E {}".format(blacklist_bed) if blacklist_bed else "",
                control_pfx=lambda wildcards,input: "-c" if input.control else "",
                control=lambda wildcards,input: ",".join(input.control) if input.control else "",
                frag_size=fragmentLength,
                ignoreForNorm = "-e " + ','.join(ignoreForNormalization) if ignoreForNormalization else ""
            conda: CONDA_CHIPSEQ_ENV
            shell: """
                Genrich -t {params.bams} {params.control_pfx} {params.control} -o {output} -r {params.blacklist} {params.ignoreForNorm} -w {params.frag_size} 2> {log}
                """


rule prep_bedgraph:
    input: "filtered_bam/{sample}.namesorted.bam"
    output: temp("filtered_bedgraph/{sample}.fragments.bedgraph")
    log: "filtered_bedgraph/log/{sample}.log"
    params:
        sample = lambda wildcards: wildcards.sample,
        genome = genome_index
    conda: CONDA_RNASEQ_ENV
    shell: """
        bedtools bamtobed -bedpe -i {input} | awk '$1==$4 && $6-$2 < 1000 {{print $0}}' - | cut -f 1,2,6 - | sort -k1,1 -k2,2n -k3,3n > filtered_bedgraph/{params.sample}.fragments.bed
        bedtools genomecov -bg -i filtered_bedgraph/{params.sample}.fragments.bed -g {params.genome} > {output}
        """

rule SEACR_peaks_stringent:
    input:
        chip = "filtered_bedgraph/{chip_sample}.fragments.bedgraph",
        control = lambda wildcards: "filtered_bedgraph/"+get_control(wildcards.chip_sample)+".fragments.bedgraph" if get_control(wildcards.chip_sample)
                 else []
    output:
        "SEACR/{chip_sample}.filtered.stringent.bed"
    log: "SEACR/logs/{chip_sample}.log"
    params:
        fdr = lambda wildcards,input: fdr if not input.control else "",
        prefix = os.path.join(outdir,"SEACR/{chip_sample}.filtered"),
        script=os.path.join(maindir, "shared","tools/SEACR-1.3/SEACR_1.3.sh")
    conda: CONDA_SEACR_ENV
    shell: """
        bash {params.script} {input.chip} {input.control} {params.fdr} "norm" "stringent" {params.prefix} 2>{log}
        """

rule SEACR_peaks_lenient:
    input:
        chip = "filtered_bedgraph/{chip_sample}.fragments.bedgraph",
        control = lambda wildcards: "filtered_bedgraph/"+get_control(wildcards.chip_sample)+".fragments.bedgraph" if get_control(wildcards.chip_sample)
                 else []
    output:
        "SEACR/{chip_sample}.filtered.lenient.bed"
    log: "SEACR/logs/{chip_sample}.log"
    params:
        fdr = lambda wildcards,input: fdr if not input.control else "",
        prefix = os.path.join(outdir,"SEACR/{chip_sample}.filtered"),
        script=os.path.join(maindir, "shared","tools/SEACR-1.3/SEACR_1.3.sh")
    conda: CONDA_SEACR_ENV
    shell: """
        bash {params.script} {input.chip} {input.control} {params.fdr} "norm" "lenient" {params.prefix} 2>{log}
        """

rule SEACR_peak_stringent_qc:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        peaks = "SEACR/{sample}.filtered.stringent.bed"
    output:
        qc = "SEACR/{sample}.filtered.stringent_peaks.qc.txt"
    params:
        genome_index = genome_index
    benchmark:
        "SEACR/.benchmark/SEACR_peak_stringent_qc.{sample}.filtered.benchmark"
    conda: CONDA_SHARED_ENV
    shell: """
        # get the number of peaks
        peak_count=`wc -l < {input.peaks}`
        
        # get the number of mapped reads
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


rule SEACR_peak_lenient_qc:
    input:
        bam = "filtered_bam/{sample}.filtered.bam",
        peaks = "SEACR/{sample}.filtered.lenient.bed"
    output:
        qc = "SEACR/{sample}.filtered.lenient_peaks.qc.txt"
    params:
        genome_index = genome_index
    benchmark:
        "SEACR/.benchmark/SEACR_peak_lenient_qc.{sample}.filtered.benchmark"
    conda: CONDA_SHARED_ENV
    shell: """
        # get the number of peaks
        peak_count=`wc -l < {input.peaks}`
        
        # get the number of mapped reads
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
