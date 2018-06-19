## get restriction site bed files
rule get_restrictionSite:
    input:
        genome_fasta
    output:
        enzyme + ".bed"
    params:
        res_seq = get_restriction_seq(enzyme)
    log:
        out = "log/get_restrictionSite.out",
        err = "log/get_restrictionSite.err"
    conda: CONDA_HIC_ENV
    shell:
        "findRestSite -f {input} --searchPattern {params.res_seq} -o {output} > {log.out} 2> {log.err}"


# Map
rule map_fastq_single_end:
    input: fastq_dir+"/{sample}{read}.fastq.gz"
    output:
        out =  "BWA/{sample}{read}.bam"
    log:
        out = "BWA/logs/{sample}{read}.out",
        err = "BWA/logs/{sample}{read}.err"
    threads: 15
    conda: CONDA_HIC_ENV
    shell:
        "echo 'mapping {input}' > {log.out} && "
        "bwa mem -A1 -B4  -E50 -L0 "
        "-t {threads} " + bwa_index + " {input}  2> {log.err} | "
        "samtools view -Shb - > {output.out}  2>> {log.err}"

## Make HiC Matrix
if(RF_resolution is True):
    rule build_matrix:
        input:
            R1 = "BWA/{sample}"+reads[0]+".bam",
            R2 = "BWA/{sample}"+reads[1]+".bam",
            bed = enzyme + ".bed"
        output:
             matrix ="HiC_matrices/{sample}_"+matrixFile_suffix+".h5",
        params:
             QCfolder="HiC_matrices/QCplots/{sample}_QC/",
             res_seq = get_restriction_seq(enzyme),
             dang_seq = get_dangling_seq(enzyme),
             region = lambda wildcards: "--region " + restrict_region if restrict_region else "",
             min_dist = MIN_RS_DISTANCE,
             max_dist = MAX_RS_DISTANCE
        log:
            out = "HiC_matrices/logs/{sample}"+matrixFile_suffix+".out",
            err = "HiC_matrices/logs/{sample}"+matrixFile_suffix+".err"
        threads: 15
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
            "-o {output.matrix} > {log.out} 2> {log.err}"
else:
    rule build_matrix:
        input:
            R1 = "BWA/{sample}"+reads[0]+".bam",
            R2 = "BWA/{sample}"+reads[1]+".bam"
        output:
            matrix = "HiC_matrices/{sample}_"+matrixFile_suffix+".h5",
            qc = "HiC_matrices/QCplots/{sample}_QC/QC_table.txt"

        params:
            QCfolder="HiC_matrices/QCplots/{sample}_QC/",
            bin_size = bin_size,
            region = lambda wildcards: "--region " + restrict_region if restrict_region else "",
            min_dist = MIN_RS_DISTANCE,
            max_dist = MAX_RS_DISTANCE
        log:
            out = "HiC_matrices/logs/{sample}"+matrixFile_suffix+".out",
            err = "HiC_matrices/logs/{sample}"+matrixFile_suffix+".err"
        threads: 15
        conda: CONDA_HIC_ENV
        shell:
            "hicBuildMatrix -s {input.R1} {input.R2} "
            "-bs {params.bin_size} "
            "--minDistance {params.min_dist} "
            "--maxDistance {params.max_dist} "
            "--QCfolder {params.QCfolder} "
            "--threads {threads} "
            "{params.region} "
            "-o {output.matrix} > {log.out} 2> {log.err}"

## Merge the samples if asked
rule merge_matrices:
      input:
          expand("HiC_matrices/{sample}_"+matrixFile_suffix+".h5", sample=samples)
      output:
          matrix = "HiC_matrices/mergedSamples_"+matrixFile_suffix+".h5",
      log:
         out = "HiC_matrices/logs/mergedSamples_"+matrixFile_suffix+".out",
         err = "HiC_matrices/logs/mergedSamples_"+matrixFile_suffix+".err"
      conda: CONDA_HIC_ENV
      shell:
          "hicSumMatrices -m {input} -o {output.matrix} > {log.out} &> {log.err}"

## Merge the bins if asked
rule merge_bins:
     input:
         "HiC_matrices/{sample}_"+matrixFile_suffix+".h5"
     output:
         matrix = "HiC_matrices/{sample}_Mbins"+str(nbins_toMerge)+"_"+matrixFile_suffix+".h5"

     params:
         num_bins=nbins_toMerge
     log:
         out = "HiC_matrices/logs/{sample}_Mbins"+str(nbins_toMerge)+"_"+matrixFile_suffix+".out",
         err = "HiC_matrices/logs/{sample}_Mbins"+str(nbins_toMerge)+"_"+matrixFile_suffix+".err"
     conda: CONDA_HIC_ENV
     shell:
         "hicMergeMatrixBins -m {input} -nb {params.num_bins} -o {output.matrix} >{log.out} 2>{log.err} "

## diagnostic plots
rule diagnostic_plot:
    input:
        "HiC_matrices/{sample}_"+matrixFile_suffix+".h5"
    output:
        plot = "HiC_matrices/QCplots/{sample}_"+matrixFile_suffix+"_diagnostic_plot.pdf",
        mad = "HiC_matrices/QCplots/{sample}_"+matrixFile_suffix+"_mad_threshold.out"
    params:
        chr = chromosomes
    conda: CONDA_HIC_ENV
    shell:
       "hicCorrectMatrix diagnostic_plot -m {input} -o {output.plot} {params.chr} &> {output.mad} "


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
rule correct_matrix:
    input:
        matrix= "HiC_matrices/{sample}_"+matrixFile_suffix+".h5",
        correct = "HiC_matrices_corrected/logs/thresholds_{sample}_"+matrixFile_suffix+".out"
    output:
        "HiC_matrices_corrected/{sample}_"+matrixFile_suffix+".corrected.h5"
    params:
        chr = chromosomes
    conda: CONDA_HIC_ENV
    shell:
        "thresholds=$(cat \"{input.correct}\");"
        "hicCorrectMatrix correct --filterThreshold $thresholds"
        " {params.chr} -m {input.matrix} -o {output} >> {input.correct}"


## Call TADs
rule call_tads:
    input:
        "HiC_matrices_corrected/{sample}_"+matrixFile_suffix+".corrected.h5"
    output:
        "TADs/{sample}_"+matrixFile_suffix+"_boundaries.bed"
    params:
        prefix="TADs/{sample}_"+matrixFile_suffix,
        parameters=tadparams
    threads: 10
    log:
        out = "TADs/logs/{sample}_findTADs.out",
        err = "TADs/logs/{sample}_findTADs.err"
    conda:
        CONDA_HIC_ENV
    shell:
        "hicFindTADs -m {input} "
        "{params.parameters} "
        "--correctForMultipleTesting bonferroni "
        "-p {threads} "
        "--outPrefix {params.prefix} > {log.out} 2> {log.err}"

##compare matrices using hicPlotDistVsCounts
rule distvscounts:
   input:
       expand("HiC_matrices_corrected/{sample}_"+matrixFile_suffix+".corrected.h5",sample=samples)
   output:
        "HiC_matrices_corrected/dist_vs_counts.png"
   params:
        function_params = lambda wildcards: distVsCountParams if distVsCountParams else " ",
   log:
        out = "HiC_matrices_corrected/logs/dist_vs_counts.out",
        err = "HiC_matrices_corrected/logs/dist_vs_counts.err"

   conda:
       CONDA_HIC_ENV
   shell:
       "hicPlotDistVsCounts -m {input} -o {output} {params.function_params} > {log.out} 2> {log.err}"
