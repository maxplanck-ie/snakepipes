
from os import path

def multiqc_input_check(return_value):
    infiles = []
    indir = ""

    if fastqc:
        infiles.append( expand("FastQC/{sample}{read}_fastqc.html", sample = samples, read = reads) )
        indir +=" FastQC "
        if trim:
            infiles.append( expand("FastQC_trimmed/{sample}{read}_fastqc.html", sample = samples, read = reads) )
            indir += " FastQC_trimmed "
    if trim:
        infiles.append( expand(fastq_dir+"/{sample}{read}.fastq.gz", sample = samples, read = reads) )
        indir += fastq_dir + " "

    if mapping_prg == "Bowtie2":
        # pipeline is DNA-mapping
        infiles.append( expand("Bowtie2/{sample}.Bowtie2_summary.txt", sample = samples) +
                  expand("Picard_qc/AlignmentSummaryMetrics/{sample}.alignment_summary_metrics.txt", sample = samples) +
                  expand("Picard_qc/MarkDuplicates/{sample}.mark_duplicates_metrics.txt", sample = samples) )
        if paired:
            infiles.append( expand("Picard_qc/InsertSizeMetrics/{sample}.insert_size_metrics.txt", sample = samples) )
        indir += " Picard_qc "
        if qualimap:
            infiles.append( expand("Qualimap_qc/{sample}.filtered.bamqc_results.txt", sample = samples) )
            indir += " Qualimap_qc "
    else:
        # must be RNA-mapping, add files as per the mode

        if not "mapping-free" in mode:
            infiles.append( expand(mapping_prg+"/{sample}.bam", sample = samples) )
            indir += mapping_prg + " featureCounts "
            if "allelic-mapping" in mode:
                infiles.append( expand("featureCounts/{sample}.allelic_counts.txt", sample = samples) )
            else:
                infiles.append( expand("featureCounts/{sample}.counts.txt", sample = samples) )
        else:
            infiles.append( expand("Salmon/{sample}/quant.sf", sample = samples) )
            indir += " Salmon "

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
    log: "multiQC/multiQC.log"
    shell:
        multiqc_path+"multiqc -o multiQC {params.indirs} &> {log}"
