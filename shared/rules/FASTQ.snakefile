# # downsample FASTQ files
# if downsample:
#     rule FASTQ:
#         input:
#             indir+"/{sample}{read}.fastq.gz"
#         output:
#             fq_sampl = "FASTQ/{sample}{read}.downsampled.fastq.gz",
#             fq_symlink = "FASTQ/{sample}{read}.fastq.gz"
#         benchmark:
#             "FASTQ/.benchmark/FASTQ_downsample.{sample}{read}.benchmark"
#         threads: 4
#         shell:
#             "bash %s {threads} %i {input} {output.fq_sampl}" % (os.path.join(workflow_tools,"fastq_head.sh"), int(config["downsample"])) + " "
#             "&& ( [ -f {output.fq_symlink} ] || ln -s -r {output.fq_sampl} {output.fq_symlink} ) && touch -h {output.fq_symlink}"
# #  create symlinks to input files
# else:
#     rule FASTQ:
#         input:
#             indir+"/{sample}{read}.fastq.gz"
#         output:
#             "FASTQ/{sample}{read}.fastq.gz"
#         shell:
#             "( [ -f {output} ] || ln -s -r {input} {output} ) && touch -h {output}"
# downsample FASTQ files
if downsample:
    if paired:
        rule FASTQ:
            input:
                r1 = indir+"/{sample}"+reads[0]+".fastq.gz",
                r2 = indir+"/{sample}"+reads[1]+".fastq.gz"
            output:
                r1_sampl = "FASTQ/{sample}"+reads[0]+".downsampled.fastq.gz",
                r2_sampl = "FASTQ/{sample}"+reads[1]+".downsampled.fastq.gz",
                r1_symlink = "FASTQ/{sample}"+reads[0]+".fastq.gz",
                r2_symlink = "FASTQ/{sample}"+reads[1]+".fastq.gz"
            params:
                num_reads = downsample
            benchmark:
                "FASTQ/.benchmark/FASTQ_downsample.{sample}.benchmark"
            threads: 10
            shell:
                os.path.join(workflow_tools,"downsample_se_pe.sh")+" {params.num_reads} {threads} {input.r1} {output.r1_sampl} {input.r2} {output.r2_sampl}"
                " && ( ( [ -f {output.r1_symlink} ] && [ -f {output.r2_symlink} ] ) || ( ln -s -r {output.r1_sampl} {output.r1_symlink} && ln -s -r {output.r2_sampl} {output.r2_symlink} ) ) && touch -h {output.r1_symlink} {output.r2_symlink} "
    else:
        rule FASTQ:
            input:
                indir+"/{sample}.fastq.gz"
            output:
                fq_sampl = "FASTQ/{sample}.downsampled.fastq.gz",
                fq_symlink = "FASTQ/{sample}.fastq.gz"
            threads: 12
            params:
                num_reads = downsample
            shell:
                os.path.join(workflow_tools,"downsample_se_pe.sh")+" {params.num_reads} {threads} {input} {output.fq_sampl} "
                " && ( [ -f {output} ] || ( ln -s -r {output.fq_sampl} {output.fq_symlink} ) ) && touch -h {output.fq_symlink}"
else:
    rule FASTQ:
        input:
            indir+"/{sample}{read}.fastq.gz"
        output:
            "FASTQ/{sample}{read}.fastq.gz"
        shell:
            "( [ -f {output} ] || ln -s -r {input} {output} ) && touch -h {output}"

#
#            "( [ -f {output} ] || ln -s -r {input} {output} ) && touch -h {output}"
