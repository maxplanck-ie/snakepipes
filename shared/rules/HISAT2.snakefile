### HISAT2 #####################################################################

if paired:
    rule HISAT2:
        input:
            r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
            r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz"
        output:
            align_summary = "HISAT2/{sample}.HISAT2_summary.txt",
            bam = "HISAT2/{sample}.bam",
            splice = "HISAT2/{sample}.splice_sites.txt",
            met = "HISAT2/{sample}.metrics.txt",
            unconc = "HISAT2/{sample}.un-conc.fastq.gz",
            alconc = "HISAT2/{sample}.al-conc.fastq.gz"
        params:
            hisat_opts = "",
            library_type = library_type
        benchmark:
            "HISAT2/.benchmark/HISAT2.{sample}.benchmark"
        threads: 8
        shell:
            hisat2_path+"hisat2 "
            "-p {threads} "
            "{params.hisat_opts} "
            "--rna-strandness {params.library_type} "
            "-x "+hisat2_index+" "
            "-1 {input.r1} -2 {input.r2} "
            "--novel-splicesite-outfile {output.splice} "
            "--met-file {output.met} "
            "--un-conc-gz {output.unconc} "
            "--al-conc-gz {output.alconc} "
            "2> {output.align_summary} | "
            ""+samtools_path+"samtools view -Sb - | "
            ""+samtools_path+"samtools sort -m 2G -T ${{TMPDIR}}{wildcards.sample} -@ 2 -O bam - > {output.bam} "
            "&& touch {output.unconc} {output.alconc} "
else:
    rule HISAT2:
        input:
            fastq_dir+"/{sample}.fastq.gz"
        output:
            align_summary = "HISAT2/{sample}.HISAT2_summary.txt",
            bam = "HISAT2/{sample}.bam",
            splice = "HISAT2/{sample}.splice_sites.txt",
            met = "HISAT2/{sample}.metrics.txt",
            un = "HISAT2/{sample}.un.fastq.gz",
            al = "HISAT2/{sample}.al.fastq.gz"
        params:
            hisat_opts = "",
            rna_strandness = rna_strandness
        benchmark:
            "HISAT2/.benchmark/HISAT2.{sample}.benchmark"
        threads: 8
        shell:
            hisat2_path+"hisat2 "
            "-p {threads} "
            "{params.hisat_opts} "
            "{params.rna_strandness} "
            "-x "+hisat2_index+" "
            "-U {input} "
            "--novel-splicesite-outfile {output.splice} "
            "--met-file {output.met} "
            "--un-gz {output.un} "
            "--al-gz {output.al} "
            "2> {output.align_summary} | "
            ""+samtools_path+"samtools view -Sb - | "
            ""+samtools_path+"samtools sort -m 2G -T ${{TMPDIR}}{wildcards.sample} -@ 2 -O bam - > {output.bam} "
            "&& touch {output.un} {output.al} "


### samtools_index #############################################################

rule HISAT2_BAM_index:
    input:
        "HISAT2/{sample}.bam"
    output:
        "HISAT2/{sample}.bam.bai"
    shell:
        samtools_path+"samtools index {input}"
