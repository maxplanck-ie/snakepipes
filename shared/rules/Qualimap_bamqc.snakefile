### Qualimap bamqc #############################################################

rule Qualimap_bamqc:
    input:
        "filtered_bam/{sample}.bam"
    output:
        txt = "Qualimap_qc/{sample}/genome_results.txt",
        html = "Qualimap_qc/{sample}/qualimapReport.html"
    params:
        outdir = "Qualimap_qc/{sample}",
        collect_overlap_pairs = "--collect-overlap-pairs" if paired else ""
    log:
        "Qualimap_qc/logs/Qualimap_bamqc.{sample}.log"
    benchmark:
        "Qualimap_qc/.benchmark/Qualimap_bamqc.{sample}.benchmark"
    threads: 16
    shell:
        "export PATH="+R_path+":$PATH && "
        # unset DISPLAY environment variable as Java VM might fail otherwise
        "unset DISPLAY && "
        ""+qualimap_path+"qualimap bamqc "
        "--java-mem-size=8G "
        "--bam {input} "
        "--paint-chromosome-limits "
        "{params.collect_overlap_pairs} "
        "-nt {threads} "
        "--outdir {params.outdir} "
        "&> {log}"


rule Qualimap_bamqc_symlink_txt:
    input:
        "Qualimap_qc/{sample}/genome_results.txt"
    output:
        "Qualimap_qc/{sample}.bamqc_results.txt"
    shell:
        "( [ -f {output} ] || ln -s -r {input} {output} ) && touch -h {output}"


rule Qualimap_bamqc_symlink_html:
    input:
        "Qualimap_qc/{sample}/qualimapReport.html"
    output:
        "Qualimap_qc/{sample}.bamqc_report.html"
    shell:
        "( [ -f {output} ] || ln -s -r {input} {output} ) && touch -h {output}"
