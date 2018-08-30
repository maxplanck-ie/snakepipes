def downloadFile(url, output):
    import urllib.request
    import gzip
    import bz2
    import os.path

    if os.path.exists(url):
        url = "file://{}".format(urllib.request.pathname2url(url))

    f = urllib.request.urlopen(url)
    content = f.read()
    f.close()

    of = open(output[0], "wb")

    # Sniff the file format and decompress as needed
    first3 = bytes(content[:3])
    if first3 == b"\x1f\x8b\x08":
        of.write(gzip.decompress(content))
    elif first3 == b"\x42\x5a\x68":
        of.write(bz2.decompress(content))
    else:
        of.write(content)
    of.close()


# Default memory allocation: 20G
rule createGenomeFasta:
    output: genome_fasta
    params:
        url = genomeURL
    run:
        downloadFile(params.url, output)


# Default memory allocation: 1G
rule fastaIndex:
    input: genome_fasta
    output: genome_index
    conda: CONDA_SHARED_ENV
    shell: """
        samtools faidx {input}
        """


# Default memory allocation: 8G
rule make2bit:
    input: genome_fasta
    output: genome_2bit
    conda: CONDA_CREATE_INDEX_ENV
    shell: """
        faToTwoBit {input} {output}
        """


# This is the same as createGenomeFasta, we could decrease this to an external script
# Default memory allocation: 20G
rule downloadGTF:
    output: genes_gtf
    params: 
        url = gtfURL
    run:
        downloadFile(params.url, output)


# Default memory allocation: 1G
rule gtf2BED:
    input: genes_gtf
    output: genes_bed
    conda: CONDA_CREATE_INDEX_ENV
    shell: """
        awk '{{if ($3 != "gene") print $0;}}' {input} \
            | grep -v "^#" \
            | gtfToGenePred /dev/stdin /dev/stdout \
            | genePredToBed stdin {output}
        """

# Default memory allocation: 1G
rule extendCodingRegions:
    input: genes_gtf
    output: extended_coding_regions_gtf
    shell: """
        grep -v "^#" {input} | awk 'BEGIN{{OFS="\t"}}{{if($3 == "gene" || $3 == "transcript") {{$4 -= 500; $5 += 500; if($4 < 1) {{$4 = 1}}; print}}}}' > {output}
        """


# Default memory allocation: 20G
rule bowtie2Index:
    input: genome_fasta
    output: os.path.join(outdir, "BowtieIndex/genome.rev.2.bt2")
    params:
      basedir = os.path.join(outdir, "BowtieIndex")
    conda: CONDA_DNA_MAPPING_ENV
    threads: 10
    shell: """
        ln -s {input} {params.basedir}/genome.fa
        bowtie2-build -t {threads} {params.basedir}/genome.fa {params.basedir}/genome
        """

# Default memory allocation: 20G
rule hisat2Index:
    input: genome_fasta
    output: os.path.join(outdir, "HISAT2Index/genome.6.ht2")
    params:
      basedir = os.path.join(outdir, "HISAT2Index")
    threads: 10
    conda: CONDA_RNASEQ_ENV
    shell: """
        ln -s {input} {params.basedir}/genome.fa
        hisat2-build -q -p {threads} {params.basedir}/genome.fa {params.basedir}/genome
        """


# Default memory allocation: 1G
rule makeKnownSpliceSites:
    input: genes_gtf
    output: known_splicesites
    conda: CONDA_RNASEQ_ENV
    threads: 10
    shell: """
        extract_splice_sites.py {input} > {output}
        """


# Default memory allocation: 80G
rule starIndex:
    input: genome_fasta
    output: os.path.join(outdir, "STARIndex/SAindex")
    params:
      basedir = os.path.join(outdir, "STARIndex")
    conda: CONDA_RNASEQ_ENV
    threads: 10
    shell: """
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.basedir} --genomeFastaFiles {input}
        rm Log.out
        """


# Default memory allocation: 8G
rule bwaIndex:
    input: genome_fasta
    output: os.path.join(outdir, "BWAIndex/genome.fa.sa")
    params:
      genome = os.path.join(outdir, "BWAIndex", "genome.fa")
    conda: CONDA_HIC_ENV
    shell: """
        ln -s {input} {params.genome}
        bwa index {params.genome}
        """


# Default memory allocation: 8G
rule bwamethIndex:
    input: genome_fasta
    output: os.path.join(outdir, "BWAmethIndex/genome.fa.bwameth.c2t.sa")
    params:
      genome = os.path.join(outdir, "BWAmethIndex", "genome.fa")
    conda: CondaEnvironment
    shell: """
        ln -s {input} {params.genome}
        bwameth.py index {params.genome}
        """


# Default memory allocation: 1G
rule copyBlacklist:
    output: os.path.join(outdir, "annotation/blacklist.bed")
    params: 
        url = blacklist
    run:
        downloadFile(params.url, output)


# Default memory allocation: 1G
rule computeEffectiveGenomeSize:
    input: genome_fasta
    output: os.path.join(outdir, "genome_fasta", "effectiveSize")
    conda: CONDA_SHARED_ENV
    shell: """
        seqtk comp {input} | awk '{{tot += $3 + $4 + $5 + $6}}END{{print tot}}' > {output}
        """
