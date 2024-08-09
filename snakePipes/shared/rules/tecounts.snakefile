import os
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
        threads: lambda wildcards: 20 if 20<max_thread else max_thread  # 3.2G per core
        conda: CONDA_RNASEQ_ENV
        shell: """
            ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
            STAR --runThreadN {threads} \
                {params.alignerOptions} \
                --sjdbOverhang 100 \
                --outSAMunmapped Within \
                --outSAMtype BAM Unsorted \
                --sjdbGTFfile {params.gtf} \
                --genomeDir {params.index} \
                --readFilesIn <(gunzip -c {input.r1}) <(gunzip -c {input.r2}) \
                --outFileNamePrefix {params.prefix}
            mv {params.prefix}Aligned.out.bam {output.bam}
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
        threads: lambda wildcards: 20 if 20<max_thread else max_thread  # 3.2G per core
        conda: CONDA_RNASEQ_ENV
        shell: """
            ( [ -d {params.sample_dir} ] || mkdir -p {params.sample_dir} )
            STAR --runThreadN {threads} \
                {params.alignerOptions} \
                --sjdbOverhang 100 \
                --outSAMunmapped Within \
                --outSAMtype BAM Unsorted \
                --sjdbGTFfile {params.gtf} \
                --genomeDir {params.index} \
                --readFilesIn <(gunzip -c {input}) \
                --outFileNamePrefix {params.prefix}
            mv {params.prefix}Aligned.out.bam {output.bam}
            """


rule makeRMSKGTF:
    input: rmsk_file
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

            of.write("{genoName}\trmsk\texon\t{genoStart}\t{genoEnd}\t.\t{strand}\t.\tgene_id \"{repName}\"; transcript_id \"{tid}\"; family_id \"{repFamily}\"; class_id \"{repClass}\";\n".format(**meta))
        f.close()
        of.close()


# Run TE Counts single core, 30GB
rule TEcounts:
    input:
        bam = aligner + "/{sample}.unsorted.bam",
        repeatGTF = "rmsk.gtf"
    output:
        "TEcount/{sample}.cntTable"
    params:
        gtf = genes_gtf
    benchmark:
        "TEcount/.benchmark/{sample}.benchmark"
    threads: 1
    conda: CONDA_NONCODING_RNASEQ_ENV
    shell: """
        TEcount --format BAM --mode multi -b {input.bam} --GTF {params.gtf} --TE {input.repeatGTF} --project TEcount/{wildcards.sample}
        """


rule sortBams:
    input:
        aligner + "/{sample}.unsorted.bam"
    output:
        "filtered_bam/{sample}.filtered.bam"
    threads: 5
    params:
        tempDir = tempDir
    conda: CONDA_SHARED_ENV
    shell: """
        TMPDIR={params.tempDir}
        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
        samtools view -u -F 2304 {input} | samtools sort -@ 4 -m 2G -T $MYTEMP/{wildcards.sample} -o {output}
        rm -rf $MYTEMP
        """
if fromBAM:
    rule samtools_index_filtered_bam:
            input:
                "filtered_bam/{sample}.filtered.bam"
            output:
                "filtered_bam/{sample}.filtered.bam.bai"
            conda: CONDA_SHARED_ENV
            shell: "samtools index {input}"


rule cpGTF:
    input:
        genes_gtf,
        genes_bed
    output:
        temp("Annotation/genes.filtered.gtf"),
        temp("Annotation/genes.filtered.bed")
    run:
        if not os.path.exists(os.path.join(outdir,output[0])):
            os.symlink(input[0],os.path.join(outdir,output[0]))
        if not os.path.exists(os.path.join(outdir,output[1])):
            os.symlink(input[1],os.path.join(outdir,output[1]))


rule symbolFile:
    input: genes_gtf
    output: "Annotation/genes.filtered.symbol"
    run:
        of = open(output[0], "w")
        for line in open(input[0]):
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            geneID = None
            geneSymbol = None
            if cols[2] == "gene":
                meta = cols[8].split(";")
                meta = [x.strip() for x in meta]  # remove random spaces
                for x in meta:
                    if x.startswith("gene_id "):
                        geneID = x[8:].strip('"')
                    elif x.startswith("gene_name "):
                        geneSymbol = x[10:].strip('"')
                if geneID is not None and geneSymbol is not None:
                    of.write("{}\t{}\n".format(geneID, geneSymbol))
        of.close()


def get_outdir(folder_name,sampleSheet):
    sample_name = os.path.splitext(os.path.basename(str(sampleSheet)))[0]
    return("{}_{}".format(folder_name, sample_name))

if sampleSheet:
    rule split_sampleSheet:
        input:
            sampleSheet = sampleSheet
        output:
            splitSheets = os.path.join("splitSampleSheets",os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv")
        params:
            splitSheetPfx = os.path.join("splitSampleSheets",os.path.splitext(os.path.basename(str(sampleSheet)))[0])
        run:
            if isMultipleComparison:
                cf.splitSampleSheet(input.sampleSheet,params.splitSheetPfx)


# TODO: topN, FDR
if sampleSheet:
    if not isMultipleComparison:
        rule DESeq2:
            input:
                cnts=["TEcount/{}.cntTable".format(x) for x in samples],
                sampleSheet=sampleSheet,
                symbol_file = "Annotation/genes.filtered.symbol"
            output:
                "{}/DESeq2.session_info.txt".format(get_outdir("DESeq2",sampleSheet))
            benchmark:
                "{}/.benchmark/DESeq2.featureCounts.benchmark".format(get_outdir("DESeq2",sampleSheet))
            params:
                outdir = os.path.join(outdir, get_outdir("DESeq2", sampleSheet)),
                fdr = 0.05,
            conda: CONDA_RNASEQ_ENV
            script: "../rscripts/noncoding-DESeq2.R"
    else:
        rule DESeq2:
            input:
                cnts=["TEcount/{}.cntTable".format(x) for x in samples],
                sampleSheet = "splitSampleSheets/" + os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv",
                symbol_file = "Annotation/genes.filtered.symbol"
            output:
                "{}/DESeq2.session_info.txt".format(get_outdir("DESeq2",os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv"))
            benchmark:
                "{}/.benchmark/DESeq2.featureCounts.benchmark".format(get_outdir("DESeq2",os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv"))
            params:
                outdir = lambda wildcards,input: os.path.join(outdir,get_outdir("DESeq2",input.sampleSheet)),
                fdr = 0.05,
            conda: CONDA_RNASEQ_ENV
            script: "../rscripts/noncoding-DESeq2.R"
