
import os

## get restriction site bed files
rule get_restrictionSite:
    input:
        genome_fasta
    output:
        enzyme + ".bed"
    params:
        res_seq = get_restriction_seq(enzyme)
    shell:
        hicExplorer_path + "findRestSite -f {input} --searchPattern {params.res_seq} -o {output}"

# Map
rule map_fastq_single_end:
    input: fastq_dir+"/{sample}{read}.fastq.gz"
    output: "BWA/{sample}{read}.bam"
    log:    "BWA/{sample}{read}.log"
    threads: 16
    shell:
        "echo 'mapping {input}' > {log} && "
        + bwa_path + "bwa mem -A1 -B4  -E50 -L0 "
        "-t {threads} " + bwa_index + " {input} 2>> {log} | "
        + samtools_path + "samtools view -Shb - > {output}"

## Make HiC Matrix
if(RF_resolution is True):
    rule build_matrix:
        input:
            R1 = "BWA/{sample}"+reads[0]+".bam",
            R2 = "BWA/{sample}"+reads[1]+".bam",
            bed = enzyme + ".bed"
        output:
             matrix = "HiC_matrices/{sample}_"+matrixFile_suffix+".h5",
             bam = "BWA/{sample}_R12_"+matrixFile_suffix+".bam"
        params:
             QCfolder="HiC_matrices/QCplots/{sample}_QC/",
             res_seq = get_restriction_seq(enzyme),
             dang_seq = get_dangling_seq(enzyme),
             region = lambda wildcards: "--region " + restrict_region if restrict_region else "",
             min_dist = MIN_RS_DISTANCE,
             max_dist = MAX_RS_DISTANCE
        log:
           "HiC_matrices/logs/{sample}.log"
        threads: 15
        shell:
            hicExplorer_path + "hicBuildMatrix -s {input.R1} {input.R2} "
            "-rs {input.bed} "
            "--restrictionSequence {params.res_seq} "
            "--danglingSequence {params.dang_seq} "
            "--minDistance {params.min_dist} "
            "--maxDistance {params.max_dist} "
            "--QCfolder {params.QCfolder} "
            "--threads {threads} "
            "{params.region} "
            "-b {output.bam} -o {output.matrix} &> {log}"
else:
    rule build_matrix:
        input:
            R1 = "BWA/{sample}"+reads[0]+".bam",
            R2 = "BWA/{sample}"+reads[1]+".bam"
        output:
            matrix = "HiC_matrices/{sample}_"+matrixFile_suffix+".h5",
            bam = "BWA/{sample}_R12_"+matrixFile_suffix+".bam"
        params:
            QCfolder="HiC_matrices/QCplots/{sample}_QC/",
            bin_size = bin_size,
            region = lambda wildcards: "--region " + restrict_region if restrict_region else "",
            min_dist = MIN_RS_DISTANCE,
            max_dist = MAX_RS_DISTANCE
        log:
           "HiC_matrices/logs/{sample}.log"
        threads: 15
        shell:
            hicExplorer_path + "hicBuildMatrix -s {input.R1} {input.R2} "
            "-bs {params.bin_size} "
            "--minDistance {params.min_dist} "
            "--maxDistance {params.max_dist} "
            "--QCfolder {params.QCfolder} "
            "--threads {threads} "
            "{params.region} "
            "-b {output.bam} -o {output.matrix} &> {log}"

## Merge the samples if asked
if(merge_samples is True):
    rule merge_matrices:
        input:
            expand("HiC_matrices/{sample}_"+matrixFile_suffix+".h5", sample=samples)
        output:
            "HiC_matrices/all_matrices_merged.h5"
        shell:
            hicExplorer_path + "hicSumMatrices -m {input} -o {output}"

## Merge the bins if asked
if(nbins_toMerge != 0):
    rule merge_bins:
        input:
            "HiC_matrices/{sample}_"+matrixFile_suffix+".h5"
        output:
            "HiC_matrices/{sample}_"+matrixFile_suffix+"_m"+str(nbins_toMerge)+".h5"
        params:
            num_bins=nbins_toMerge
        shell:
            hicExplorer_path + "hicMergeMatrixBins -m {input} -nb {params.num_bins} -o {output}"

## diagnostic plots
rule diagnostic_plot:
    input:
        "HiC_matrices/{sample}_"+matrixFile_suffix+".h5"
    output:
        plot = "HiC_matrices/QCplots/{sample}_"+matrixFile_suffix+"_diagnostic_plot.pdf",
        mad = "HiC_matrices/QCplots/{sample}_"+matrixFile_suffix+"_mad_threshold.out"
    shell:
        hicExplorer_path + "hicCorrectMatrix diagnostic_plot -m {input} -o {output.plot} > {output.mad}"

## Correct matrices
rule correct_matrix:
    input:
        matrix = "HiC_matrices/{sample}_"+matrixFile_suffix+".h5",
        mad = "HiC_matrices/QCplots/{sample}_"+matrixFile_suffix+"_mad_threshold.out"
    output:
        "HiC_matrices_corrected/{sample}_"+matrixFile_suffix+".corrected.h5"
    log:
        correct = "HiC_matrices_corrected/logs/{sample}_"+matrixFile_suffix+".log"
    run:
        thresholds = get_mad_score(input.mad)
        f = open(log.correct, 'w')
        f.write('Thresholds for matrix correction are : {} \n'. format(thresholds))
        f.close()

        shell(
            hicExplorer_path + "hicCorrectMatrix correct --filterThreshold " +
            thresholds + " -m {input.matrix} -o {output} >> {log.correct} 2>&1"
            )


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
       "TADs/logs/{sample}_findTADs.log"
    shell:
        hicExplorer_path + "hicFindTADs -m {input} "
        "{params.parameters} "# needs to be variable
        "--correctForMultipleTesting bonferroni "
        "-p {threads} "
        "--outPrefix {params.prefix} > {log} 2>&1"
