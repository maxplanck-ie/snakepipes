def getAlignerCmd(which_aligner):
    cmd_str=""
    if which_aligner=="bwa":
        cmd_str="bwa mem"
    elif which_aligner=="bwa-mem2":
        cmd_str="bwa-mem2 mem"
    return(cmd_str)

def getAlignerIndex(which_aligner):
    which_index=""
    if which_aligner=="bwa":
        which_index=bwa_index
    elif which_aligner=="bwa-mem2":
        which_index=bwa_mem2_index
    return(which_index)

## get restriction site bed files
rule get_restrictionSite:
    input:
        genome_fasta
    output:
        enzyme + ".bed"
    params:
        res_seq = get_restriction_seq(enzyme)
    conda: CONDA_HIC_ENV
    shell:
        "hicFindRestSite -f {input} --searchPattern {params.res_seq} -o {output}"


# Map
rule map_fastq_single_end:
    input: fastq_dir+"/{sample}{read}.fastq.gz"
    output:
        out =  aligner+"/{sample}{read}.bam"
    params:
        aligner_cmd = getAlignerCmd(aligner),
        aligner_index = getAlignerIndex(aligner)
    threads: lambda wildcards: 15 if 15<max_thread else max_thread
    conda: CONDA_HIC_ENV
    shell:
        "echo 'mapping {input}' && "
        "{params.aligner_cmd} -A1 -B4  -E50 -L0 "
        "-t {threads} {params.aligner_index} {input} | "
        "samtools view -Shb - > {output.out}"

## Make HiC Matrix
if(RFResolution is True):
    rule build_matrix:
        input:
            R1 = aligner+"/{sample}"+reads[0]+".bam",
            R2 = aligner+"/{sample}"+reads[1]+".bam",
            bed = enzyme + ".bed"
        output:
             matrix ="HiC_matrices/{sample}_"+matrixFile_suffix+matrix_format,
             qc = "HiC_matrices/QCplots/{sample}_QC/QC.log"
        params:
             QCfolder="HiC_matrices/QCplots/{sample}_QC/",
             res_seq = get_restriction_seq(enzyme),
             dang_seq = get_dangling_seq(enzyme),
             region = lambda wildcards: "--region " + str(restrictRegion) if restrictRegion else "",
             min_dist = MIN_RS_DISTANCE,
             max_dist = MAX_RS_DISTANCE
        threads: lambda wildcards: 10 if 10<max_thread else max_thread
        conda: CONDA_HIC_ENV
        shell:
            "hicBuildMatrix -s {input.R1} {input.R2} "
            "-rs {input.bed} "
            "--restrictionSequence {params.res_seq} "
            "--danglingSequence {params.dang_seq} "
            "--minDistance {params.min_dist} "
            "--maxDistance {params.max_dist} "
            "--QCfolder {params.QCfolder} "
            "--threads {threads} "
            "{params.region} "
            "-o {output.matrix} &&"
            " rm {params.QCfolder}"+"QC_table.txt"
else:
    rule build_matrix:
        input:
            R1 = aligner+"/{sample}"+reads[0]+".bam",
            R2 = aligner+"/{sample}"+reads[1]+".bam",
            bed = enzyme + ".bed"
        output:
            matrix = "HiC_matrices/{sample}_"+matrixFile_suffix+matrix_format,
            qc = "HiC_matrices/QCplots/{sample}_QC/QC.log"

        params:
            QCfolder="HiC_matrices/QCplots/{sample}_QC/",
            bin_size = binSize,
            res_seq = get_restriction_seq(enzyme),
            dang_seq = get_dangling_seq(enzyme),
            region = lambda wildcards: "--region " + str(restrictRegion) if restrictRegion else "",
            min_dist = MIN_RS_DISTANCE,
            max_dist = MAX_RS_DISTANCE
        threads: lambda wildcards: 10 if 10<max_thread else max_thread
        conda: CONDA_HIC_ENV
        shell:
            "hicBuildMatrix -s {input.R1} {input.R2} "
            "-bs {params.bin_size} "
            "-rs {input.bed} "
            "--restrictionSequence {params.res_seq} "
            "--danglingSequence {params.dang_seq} "
            "--minDistance {params.min_dist} "
            "--maxDistance {params.max_dist} "
            "--QCfolder {params.QCfolder} "
            "--threads {threads} "
            "{params.region} "
            "-o {output.matrix} &&"
            " rm {params.QCfolder}"+"QC_table.txt"

## Merge the samples if asked
rule merge_matrices:
      input:
          lambda wildcards: expand("HiC_matrices/{sample}_"+matrixFile_suffix+matrix_format, sample = sample_dict[wildcards.group])
      output:
          matrix = "HiC_matrices/mergedSamples_{group}_"+matrixFile_suffix+matrix_format
      conda: CONDA_HIC_ENV
      shell:
          "hicSumMatrices -m {input} -o {output.matrix}"

## Merge the bins if asked
rule merge_bins:
     input:
         "HiC_matrices/{sample}_"+matrixFile_suffix+matrix_format
     output:
         matrix = "HiC_matrices/{sample}_Mbins" + str(nBinsToMerge) + "_" + matrixFile_suffix+matrix_format
     params:
         num_bins=nBinsToMerge
     conda: CONDA_HIC_ENV
     shell:
         "hicMergeMatrixBins -m {input} -nb {params.num_bins} -o {output.matrix} "

## diagnostic plots
rule diagnostic_plot:
    input:
        "HiC_matrices/{sample}_"+matrixFile_suffix+matrix_format
    output:
        plot = "HiC_matrices/QCplots/{sample}_"+matrixFile_suffix+"_diagnostic_plot.pdf",
        mad = "HiC_matrices/QCplots/{sample}_"+matrixFile_suffix+"_mad_threshold.out"
    params:
        chr = lambda wildcards: " --chromosomes " + chromosomes if chromosomes else ""
    conda: CONDA_HIC_ENV
    shell:
       "hicCorrectMatrix diagnostic_plot -m {input} -o {output.plot} {params.chr} 2> {output.mad} "


# Compute MAD score thresholds
rule compute_thresholds:
   input:
      "HiC_matrices/QCplots/{sample}_"+matrixFile_suffix+"_mad_threshold.out"
   output:
      "HiC_matrices_corrected/logs/thresholds_{sample}_"+matrixFile_suffix+".out"
   shell:
         "madscore=$(grep \"mad threshold \" {input} | sed 's/INFO:hicexplorer.hicCorrectMatrix:mad threshold //g');"
         "upper=$(echo -3*$madscore | bc);"
         "echo $madscore \" \" $upper >> {output}"


## Correct matrices
if correctionMethod == 'ICE':
    rule correct_matrix:
        input:
            matrix= "HiC_matrices/{sample}_"+matrixFile_suffix+matrix_format,
            correct = "HiC_matrices_corrected/logs/thresholds_{sample}_"+matrixFile_suffix+".out"
        output:
            "HiC_matrices_corrected/{sample}_"+matrixFile_suffix+".corrected"+matrix_format
        params:
            chr = lambda wildcards: " --chromosomes " + chromosomes if chromosomes else ""
        conda: CONDA_HIC_ENV
        shell:
            "thresholds=$(cat \"{input.correct}\");"
            "hicCorrectMatrix correct --correctionMethod ICE --filterThreshold $thresholds"
            " {params.chr} -m {input.matrix} -o {output} >> {input.correct}"

else:
     rule correct_matrix:
         input:
             matrix = "HiC_matrices/{sample}_"+matrixFile_suffix+matrix_format
         output:
             "HiC_matrices_corrected/{sample}_"+matrixFile_suffix+".corrected"+matrix_format
         params:
             chr = lambda wildcards: " --chromosomes " + chromosomes if chromosomes else ""
         conda: CONDA_HIC_ENV
         shell:
             "hicCorrectMatrix correct --correctionMethod KR "
             " {params.chr} -m {input.matrix} -o {output}"

## Call TADs
rule call_tads:
    input:
        "HiC_matrices_corrected/{sample}_"+matrixFile_suffix+".corrected"+matrix_format
    output:
        "TADs/{sample}_"+matrixFile_suffix+"_boundaries.bed"
    params:
        prefix="TADs/{sample}_"+matrixFile_suffix,
        parameters=lambda wildcards: findTADParams if findTADParams else ""
    threads: lambda wildcards: 10 if 10<max_thread else max_thread
    conda:
        CONDA_HIC_ENV
    shell:
        "hicFindTADs -m {input} "
        "{params.parameters} "
        "--correctForMultipleTesting bonferroni "
        "-p {threads} "
        "--outPrefix {params.prefix}"

##compare matrices using hicPlotDistVsCounts
rule distvscounts:
   input:
        matrices = expand("HiC_matrices_corrected/{sample}_"+matrixFile_suffix+".corrected"+matrix_format, sample = samples)
   output:
        "dist_vs_counts.png"
   params:
        function_params = lambda wildcards: distVsCountParams if distVsCountParams else " "
   conda:
       CONDA_HIC_ENV
   shell:
       "hicPlotDistVsCounts -m  {input.matrices} -o {output} {params.function_params}"
