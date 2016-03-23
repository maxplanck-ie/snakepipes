## Qualimap: bamqc #############################################################

if paired:
    rule Qualimap_bamqc:
        input:
            bam = "Bowtie2/{sample}.bam"
        output:
            txt = "Qualimap/{sample}/genome_results.txt",
            html = "Qualimap/{sample}/qualimapReport.html"
        params:
            outdir = "Qualimap/{sample}",
            collect_overlap_pairs = "--collect-overlap-pairs"
        log:    "Qualimap/log/{sample}.bamqc.log"
        benchmark:  "Qualimap/.benchmark/Qualimap_bamqc.{sample}.benchmark"
        threads: 10
        shell:  "export PATH="+R_dir+":$PATH && "
                ""+qualimap_path+" bamqc --java-mem-size=4G -nt {threads} "
                "{params.collect_overlap_pairs} --paint-chromosome-limits "
                "--bam {input.bam} --outdir {params.outdir} "
                "2>&1 | tee {log} "
else:
    rule Qualimap_bamqc:
        input:
            bam = "Bowtie2/{sample}.bam"
        output:
            txt = "Qualimap/{sample}/genome_results.txt",
            html = "Qualimap/{sample}/qualimapReport.html"
        params:
            outdir = "Qualimap/{sample}"
        log:    "Qualimap/log/{sample}.bamqc.log"
        benchmark:  "Qualimap/.benchmark/Qualimap_bamqc.{sample}.benchmark"
        threads: 10
        shell:  "export PATH="+R_dir+":$PATH && "
                ""+qualimap_path+" bamqc --java-mem-size=4G -nt {threads} "
                "--paint-chromosome-limits "
                "--bam {input.bam} --outdir {params.outdir} "
                "2>&1 | tee {log} "


rule Qualimap_bamqc_symlinks_txt:
    input:  "Qualimap/{sample}/genome_results.txt"
    output: "Qualimap/{sample}.bamqc.txt"
    shell:  "( [ -f {output} ] || ln -s -r {input} {output} ) && touch -h {output}"


rule Qualimap_bamqc_symlinks_html:
    input:  "Qualimap/{sample}/qualimapReport.html"
    output: "Qualimap/{sample}.bamqc.html"
    shell:  "( [ -f {output} ] || ln -s -r {input} {output} ) && touch -h {output}"
