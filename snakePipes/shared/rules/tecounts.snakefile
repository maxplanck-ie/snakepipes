if pairedEnd:
    rule STAR:
        input:
            r1 = fastq_dir+"/{sample}"+reads[0]+".fastq.gz",
            r2 = fastq_dir+"/{sample}"+reads[1]+".fastq.gz"
        output:
            bam = temp(aligner+"/{sample}.unsorted.bam")
        params:
            alignerOptions = str(alignerOptions or ''),
            gtf = genes_gtf,
            index = star_index,
            prefix = aligner+"/{sample}/{sample}.",
            sample_dir = aligner+"/{sample}"
        benchmark:
            aligner+"/.benchmark/STAR.{sample}.benchmark"
        threads: 20  # 3.2G per core
        conda: CONDA_RNASEQ_ENV
        shell: """
            ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
            STAR --runThreadN {threads} \
                {params.alignerOptions} \
                --sjdbOverhang 100 \
                --readFilesCommand zcat \
                --outSAMunmapped Within \
                --outSAMtype BAM Unsorted \
                --sjdbGTFfile {params.gtf} \
                --genomeDir {params.index} \
                --readFilesIn {input.r1} {input.r2} \
                --outFileNamePrefix {params.prefix} \
            mv {params.prefix}/{wildcards.sample}.bam {output.bam}
            """
else:
    rule STAR:
        input:
            fastq_dir+"/{sample}"+reads[0]+".fastq.gz"
        output:
            bam = temp(aligner+"/{sample}.unsorted.bam")
        params:
            alignerOptions = str(alignerOptions or ''),
            gtf = genes_gtf,
            index = star_index,
            prefix = aligner+"/{sample}/{sample}.",
            sample_dir = aligner+"/{sample}"
        benchmark:
            aligner+"/.benchmark/STAR.{sample}.benchmark"
        threads: 20  # 3.2G per core
        conda: CONDA_RNASEQ_ENV
        shell: """
            ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
            STAR --runThreadN {threads} \
                {params.alignerOptions} \
                --sjdbOverhang 100 \
                --readFilesCommand zcat \
                --outSAMunmapped Within \
                --outSAMtype BAM Unsorted \
                --sjdbGTFfile {params.gtf} \
                --genomeDir {params.index} \
                --readFilesIn {input} \
                --outFileNamePrefix {params.prefix} \
            mv {params.prefix}/{wildcards.sample}.bam {output.bam}
            """


rule makeRMSKGTF:
    input: TEGTF
    output: temp("rmsk.gtf")
    run:
        f = open(input[0])
        of = open(output[0], "w")
        found = dict()

        for line in f:
            if line.startswith("#"):
                continue
            bin, swScore, milliDiv, milliDel, milliIns,  genoName, genoStart, genoEnd, genoLeft, strand, repName, repClass, repFamily, repStart, repEnd, repLeft, id = line.strip().split("\t")
            if repName not in found:
                found[repName] = 0
                tid = repName
            else:
                found[repName] += 1
                tid = "{}_dup{}".format(repName, found[repName])

            meta = {"genoName": genoName,
                    "genoStart": genoStart,
                    "genoEnd": genoEnd,
                    "strand": strand,
                    "repName": repName,
                    "repClass": repClass,
                    "repFamily": repFamily,
                    "tid": tid}

            of.write("{genoName}\trmsk\texon\t{genoStart}\t{genoEnd}\t.\{strand}\t.\tgene_id \"{repName}\"; transcript_id \"{tid}\"; family_id \"{repFamily}\"; class_id \"{repClass}\";\".format(**meta))
        f.close()
        of.close()


# Run TE Counts single core, 30GB
rule TEcounts:
    input:
        bam = aligner + "/{sample}.unsorted.bam",
        repeatGTF = "rmsk.gtf"
    output:
        "TEcount/{sample}/counts.txt"
    params:
        gtf = genes_gtf
    log:
        out = "TEcount/logs/{sample}.out",
        err = "TEcount/logs/{sample}.err"
    benchmark:
        "TEcount/.benchmark/{sample}.benchmark"
    threads: 1
    conda: CONDA_RIBOMINUS_RNASEQ_ENV
    shell: """
        TEcount --format BAM --mode multi -b {input.bam} --GTF {params.gtf} --TE {input.repeatGTF} --project TEcount/{wildcards.sample} 2> {log.err} > {log.out}
        """
