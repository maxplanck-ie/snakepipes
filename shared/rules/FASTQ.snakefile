# downsample FASTQ files
if downsample:
    rule FASTQ:
        input:
            indir+"/{sample}{read}.fastq.gz"
        output:
            fq_sampl = "FASTQ/{sample}{read}.downsampled.fastq.gz",
            fq_symlink = "FASTQ/{sample}{read}.fastq.gz"
        benchmark:
            "FASTQ/.benchmark/FASTQ_downsample.{sample}{read}.benchmark"
        threads: 4
        shell:
            "bash %s {threads} %i {input} {output.fq_sampl}" % (os.path.join(workflow_tools,"fastq_head.sh"), int(config["downsample"])) +""
            "&& ( [ -f {output.fq_symlink} ] || ln -s -r {output.fq_sampl} {output.fq_symlink} ) && touch -h {output.fq_symlink}"
#  create symlinks to input files
else:
    rule FASTQ:
        input:
            indir+"/{sample}{read}.fastq.gz"
        output:
            "FASTQ/{sample}{read}.fastq.gz"
        shell:
            "( [ -f {output} ] || ln -s -r {input} {output} ) && touch -h {output}"
