### Qualimap bamqc #############################################################

rule Qualimap_bamqc:
    input:
        "filtered_bam/{sample}.filtered.bam"
    output:
        txt = "Qualimap_qc/{sample}.filtered/genome_results.txt",
        html = "Qualimap_qc/{sample}.filtered/qualimapReport.html"
    params:
        outdir = "Qualimap_qc/{sample}.filtered",
        collect_overlap_pairs = "--collect-overlap-pairs" if pairedEnd else ""
    benchmark:
        "Qualimap_qc/.benchmark/Qualimap_bamqc.{sample}.filtered.benchmark"
    threads: lambda wildcards: 16 if 16<max_thread else max_thread
    conda: CONDA_DNA_MAPPING_ENV
    shell:
        "unset DISPLAY && "
        "qualimap bamqc "
        "--java-mem-size=8G "
        "--bam {input} "
        "--paint-chromosome-limits "
        "{params.collect_overlap_pairs} "
        "-nt {threads} "
        "--outdir {params.outdir} "


rule Qualimap_bamqc_symlink_txt:
    input:
        "Qualimap_qc/{sample}.filtered/genome_results.txt"
    output:
        "Qualimap_qc/{sample}.filtered.bamqc_results.txt"
    shell: """
         ln -s ../{input} {output}
            """


rule Qualimap_bamqc_symlink_html:
    input:
        "Qualimap_qc/{sample}.filtered/qualimapReport.html"
    output:
        "Qualimap_qc/{sample}.filtered.bamqc_report.html"
    shell: """
         ln -s ../{input} {output}
            """
