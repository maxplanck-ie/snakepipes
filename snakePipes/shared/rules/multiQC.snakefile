
from os import path

def multiqc_input_check(return_value):
    infiles = []
    indir = ""
    readsIdx = 1
    if pairedEnd:
        readsIdx = 2

    if not pipeline=="scrnaseq" and ("fromBAM" not in globals() or not fromBAM):
        if pairedEnd:
            if trim and fastqc:
                infiles.append( expand("FastQC_trimmed/{sample}{read}_fastqc.html", sample = samples, read = reads) )
                indir += " FastQC_trimmed "
                infiles.append( expand(fastq_dir+"/{sample}{read}.fastq.gz", sample = samples, read = reads) )
                indir += fastq_dir + " "
            elif fastqc:
                infiles.append( expand("FastQC/{sample}{read}_fastqc.html", sample = samples, read = reads) )
                indir +=" FastQC "
        else:
            if trim and fastqc:
                infiles.append( expand("FastQC_trimmed/{sample}"+reads[0]+"_fastqc.html", sample = samples) )
                indir += " FastQC_trimmed "
                infiles.append( expand(fastq_dir+"/{sample}"+reads[0]+".fastq.gz", sample = samples) )
                indir += fastq_dir + " "
            elif fastqc:
                infiles.append( expand("FastQC/{sample}"+reads[0]+"_fastqc.html", sample = samples) )
                indir +=" FastQC "
    if pipeline=="dnamapping":
        # pipeline is DNAmapping
        if aligner=="Bowtie2":
            infiles.append("deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv")
            infiles.append(expand("Bowtie2/{sample}.Bowtie2_summary.txt", sample = samples) +
                    expand("Sambamba/{sample}.markdup.txt", sample = samples) +
                    expand("deepTools_qc/estimateReadFiltering/{sample}_filtering_estimation.txt",sample=samples))
            indir += " Sambamba "
            indir += " Bowtie2 "
            indir += " deepTools_qc "
            if qualimap:
                infiles.append( expand("Qualimap_qc/{sample}.filtered.bamqc_results.txt", sample = samples) )
                indir += " Qualimap_qc "
        elif aligner=="bwa":
            infiles.append( expand("bwa/{sample}.bwa_summary.txt", sample = samples) +
                            expand("Sambamba/{sample}.markdup.txt", sample = samples) +
                            expand("deepTools_qc/estimateReadFiltering/{sample}_filtering_estimation.txt",sample=samples))
            indir += " Sambamba "
            indir += " bwa "
            indir += " deepTools_qc "
            if qualimap:
                infiles.append( expand("Qualimap_qc/{sample}.filtered.bamqc_results.txt", sample = samples) )
                indir += " Qualimap_qc "
        elif aligner=="bwa-mem2":
            infiles.append( expand("bwa-mem2/{sample}.bwa-mem2_summary.txt", sample = samples) +
                            expand("Sambamba/{sample}.markdup.txt", sample = samples) +
                            expand("deepTools_qc/estimateReadFiltering/{sample}_filtering_estimation.txt",sample=samples))
            indir += " Sambamba "
            indir += " bwa-mem2 "
            indir += " deepTools_qc "
            if qualimap:
                infiles.append( expand("Qualimap_qc/{sample}.filtered.bamqc_results.txt", sample = samples) )
                indir += " Qualimap_qc "
        if "allelic-mapping" in mode:
            infiles.append( expand("allelic_bams/{sample}.filtered.SNPsplit_report.yaml", sample = samples) )
            infiles.append( expand("allelic_bams/{sample}.filtered.SNPsplit_sort.yaml", sample = samples) )
            indir += "allelic_bams"
    elif pipeline=="rna-seq":
        # must be RNA-mapping, add files as per the mode
        if ( "alignment" in mode or "deepTools_qc" in mode or "three-prime-seq" in mode ) and not "allelic-mapping" in mode and not "allelic-counting" in mode:
            infiles.append( expand(aligner+"/{sample}.bam", sample = samples) +
                    expand("Sambamba/{sample}.markdup.txt", sample = samples) +
                    expand("deepTools_qc/estimateReadFiltering/{sample}_filtering_estimation.txt",sample=samples)+
                    expand("featureCounts/{sample}.counts.txt", sample = samples))
            indir += aligner + " featureCounts "
            indir += " Sambamba "
            indir += " deepTools_qc "
        if "allelic-mapping" in mode or "allelic-counting" in mode:
            infiles.append( expand("featureCounts/{sample}.allelic_counts.txt", sample = samples) )
            indir += aligner + " featureCounts "
        if "allelic-mapping" in mode:
            infiles.append( expand("allelic_bams/{sample}.SNPsplit_report.yaml", sample = samples) )
            infiles.append( expand("allelic_bams/{sample}.SNPsplit_sort.yaml", sample = samples) )
            indir += "allelic_bams"
        if "alignment-free" in mode:
            if "allelic-mapping" in mode:
                infiles.append( expand("SalmonAllelic/{sample}.{allelic_suffix}/quant.sf", sample = samples,allelic_suffix=allelic_suffix) )
                indir += " SalmonAllelic "
            else:
                infiles.append( expand("Salmon/{sample}/quant.sf", sample = samples) )
                indir += " Salmon "
    elif pipeline == "ncRNAseq":
        infiles.append(expand("deepTools_qc/estimateReadFiltering/{sample}_filtering_estimation.txt",sample=samples))
        indir += " STAR deepTools_qc "
    elif pipeline == "hic":
        infiles.append(expand("HiC_matrices/QCplots/{sample}_QC/QC.log", sample = samples))
        indir += " " + aligner + " " 
        indir += " ".join(expand("HiC_matrices/QCplots/{sample}_QC ", sample = samples))
    elif pipeline == "scrnaseq":
        if trim:
            infiles.append( expand("FastQC_trimmed/{sample}"+reads[0]+"_fastqc.html", sample = samples) )
            indir += " FastQC_trimmed "
        else:
            infiles.append( expand("FastQC/{sample}{read}_fastqc.html", sample = samples, read = reads) )
            indir +=" FastQC "
        if mode == "Gruen":
            infiles.append( expand(fastq_dir+"/{sample}"+reads[0]+".fastq.gz", sample = samples) +
                            expand("Counts/{sample}.summary", sample = samples) )
            indir += fastq_dir + " Counts "
            infiles.append( expand(aligner+"/{sample}.bam", sample = samples) +
            expand("Sambamba/{sample}.markdup.txt", sample = samples) +
            expand("deepTools_qc/estimateReadFiltering/{sample}_filtering_estimation.txt", sample=samples))
            indir += aligner
            indir += " Sambamba "
            indir += " deepTools_qc "
        elif mode == "STARsolo":
            infiles.append( expand(fastq_dir+"/{sample}"+reads[0]+".fastq.gz", sample = samples) )
            infiles.append( expand(aligner+"/{sample}.bam", sample = samples) +
            expand("Sambamba/{sample}.markdup.txt", sample = samples) +
            expand("deepTools_qc/estimateReadFiltering/{sample}_filtering_estimation.txt", sample=samples))
            indir += aligner
            indir += " Sambamba "
            indir += " deepTools_qc "
        elif mode == "Alevin":
            infiles.append( expand("multiQC/Alevin_{sample}.html", sample = samples))
            indir += " Alevin "
    elif pipeline == "WGBS":
        infiles.append( expand("QC_metrics/{sample}.flagstat", sample = samples) )
        indir += " QC_metrics"
    elif pipeline == "preprocessing":
        if fastqc and optDedupDist > 0:
            infiles.append("deduplicatedFASTQ/optical_dedup_mqc.json")
            indir += " deduplicatedFASTQ"

    if return_value == "infiles":
        return(infiles)
    else:
        return(indir)


rule multiQC:
    input:
        multiqc_input_check(return_value = "infiles")
    output: "multiQC/multiqc_report.html"
    params:
        indirs = multiqc_input_check(return_value = "indir")
    log:
        out = "multiQC/multiQC.out",
        err = "multiQC/multiQC.err"
    conda: CONDA_SHARED_ENV
    shell:
        "multiqc -o multiQC -f {params.indirs} > {log.out} 2> {log.err}"
