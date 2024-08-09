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
if not spikeinGenomeURL:
    rule createGenomeFasta:
        output: genome_fasta
        params:
            url = genomeURL
        run:
            downloadFile(params.url, output)

else:
    rule createHostGenomeFasta:
        output: temp(os.path.join(outdir, "genome_fasta/host.genome.fa"))
        params:
            url = genomeURL
        run:
            downloadFile(params.url, output)

    rule createSpikeinGenomeFasta:
        output: temp(os.path.join(outdir, "genome_fasta/spikein.genome.fa"))
        params:
            url = spikeinGenomeURL
        run:
            downloadFile(params.url, output)

    rule renameSpikeinChromsFasta:
        input: os.path.join(outdir, "genome_fasta/spikein.genome.fa")
        output: temp(os.path.join(outdir, "genome_fasta/spikein.genome_renamed.fa"))
        params:
            spikeinExt = spikeinExt
        shell: """
            sed -r 's/[[:space:]]+/{spikeinExt} /' {input} > {output}
        """

    rule createGenomeFasta:
        input:
            host_fasta = os.path.join(outdir,"genome_fasta/host.genome.fa"),
            spikein_fasta = os.path.join(outdir,"genome_fasta/spikein.genome_renamed.fa")
        output: genome_fasta
        shell: """
            cat {input.host_fasta} {input.spikein_fasta} > {output}
        """


# Default memory allocation: 1G
rule fastaIndex:
    input: genome_fasta
    output: genome_index
    conda: CONDA_SHARED_ENV
    shell: """
        samtools faidx {input}
        """

# Default memory allocation: 4G
rule fastaDict:
    input: genome_fasta
    output: genome_dict
    conda: CONDA_SHARED_ENV
    shell: """
        samtools dict -o {output} {input}
        """

if rmsk_file:
    rule fetchRMSK:
        output: rmsk_file
        params:
            url = rmskURL
        run:
            downloadFile(params.url, output)

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

rule downloadSpikeinGTF:
    output: temp(os.path.join(outdir, "annotation/spikein_genes_ori.gtf"))
    params:
        url = spikeinGtfURL
    run:
        downloadFile(params.url, output)

rule renameSpikeinChromsGTF:
    input: os.path.join(outdir,"annotation/spikein_genes_ori.gtf")
    output: spikein_genes_gtf
    params:
        spikeinExt = spikeinExt
    shell: """
        awk -v FS='\\t' -v OFS='\\t' '{{ if($1 !~ /^#/){{$1=$1\"{params.spikeinExt}\"; print $0 }} else{{print $0}} }}' {input} > {output}
    """


# Default memory allocation: 1G
#rule gtf2BED:
#    input: genes_gtf
#    output: genes_bed
#    conda: CONDA_CREATE_INDEX_ENV
#    shell: """
#        awk '{{if ($3 != "gene") print $0;}}' {input} \
#            | grep -v "^#" \
#            | gtfToGenePred /dev/stdin /dev/stdout \
#            | genePredToBed stdin {output}
#        """


rule gtf_to_files:
    input:
        gtf = genes_gtf
    output:
        genes_t2g,
        os.path.join(outdir, "annotation/genes.symbol"),
        genes_bed
    run:
        import shlex
        import re

        t2g = open(output[0], "w")
        symbol = open(output[1], "w")
        GTFdict = dict()

        for line in open(input.gtf):
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            annos = re.split(''';(?=(?:[^'"]|'[^']*'|"[^"]*")*$)''', cols[8])
            if cols[2] == "gene":
                # get the gene_name and gene_id values
                gene_id = None
                gene_name = None
                for anno in annos:
                    anno = shlex.split(anno.strip(), " ")
                    if len(anno) == 0:
                        continue
                    if anno[0] == "gene_id":
                        gene_id = anno[1]
                    elif anno[0] == "gene_name":
                        gene_name = anno[1]
                if gene_id:
                    symbol.write("{}\t{}\n".format(gene_id, "" if not gene_name else gene_name))
            elif cols[2] == "transcript" or 'RNA' in cols[2]:
                # get the gene_id and transcript_id values
                gene_id = None
                transcript_id = None
                gene_name = ""
                for anno in annos:
                    anno = shlex.split(anno.strip(), " ")
                    if len(anno) == 0:
                        continue
                    if anno[0] == "gene_id":
                        gene_id = anno[1]
                    elif anno[0] == "transcript_id":
                        transcript_id = anno[1]
                    elif anno[0] == "gene_name":
                        gene_name = anno[1]
                if transcript_id:
                    t2g.write("{}\t{}\t{}\n".format(transcript_id, "" if not gene_id else gene_id, gene_name))
                    # chrom, start, end, strand, exon width and exon start offset
                    GTFdict[transcript_id] = [cols[0], cols[3], cols[4], cols[6], [], []]
            elif cols[2] == "exon":
                # get the transcript_id
                transcript_id = None
                for anno in annos:
                    anno = shlex.split(anno.strip(), " ")
                    if len(anno) == 0:
                        continue
                    if anno[0] == "transcript_id":
                        transcript_id = anno[1]
                if transcript_id and transcript_id in GTFdict:
                    exonWidth = int(cols[4]) - int(cols[3]) + 1
                    exonOffset = int(cols[3]) - int(GTFdict[transcript_id][1])
                    GTFdict[transcript_id][4].append(str(exonWidth))
                    GTFdict[transcript_id][5].append(str(exonOffset))

        t2g.close()
        symbol.close()

        BED = open(output[2], "w")
        for k, v in GTFdict.items():
            # sort the starts and sizes together
            v[5] = [int(x) for x in v[5]]
            v[4] = [int(x) for x in v[4]]
            blockSizes = [str(x) for _,x in sorted(zip(v[5], v[4]))]
            blockStarts = sorted(v[5])
            blockStarts = [str(x) for x in blockStarts]
            BED.write("{}\t{}\t{}\t{}\t.\t{}\t{}\t{}\t255,0,0\t{}\t{}\t{}\n".format(v[0],  # chrom
                                                                               v[1],  # start
                                                                               v[2],  # end
                                                                               k,
                                                                               v[3],  # strand
                                                                               v[1],  # start
                                                                               v[2],  # end
                                                                               len(v[4]),  # blockCount
                                                                               ",".join(blockSizes),  # blockSizes
                                                                               ",".join(blockStarts)))  # blockStarts
        BED.close()




# Default memory allocation: 1G
# As a side effect, this checks the GTF and fasta file for chromosome name consistency (it will pass if at least 1 chromosome name is shared)
rule extendGenicRegions:
    input: genes_gtf, genome_index
    output: extended_coding_regions_gtf
    run:
        import sys
        import os

        faiChroms = set()
        for line in open(input[1]):
            cols = line.strip().split()
            faiChroms.add(cols[0])

        gtfChroms = set()
        o = open(output[0], "w")
        for line in open(input[0]):
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            gtfChroms.add(cols[0])
            if cols[2] == "gene" or cols[2] == "transcript":
                cols[3] = str(max(1, int(cols[3]) - 500))
                cols[4] = str(int(cols[4]) + 500)
            o.write("\t".join(cols))
            o.write("\n")
        o.close()

        # Ensure there is at least one shared chromosome name between the annotation and fasta file
        try:
            assert len(faiChroms.intersection(gtfChroms)) >= 1
        except:
            os.remove(output[0])
            sys.exit("There are no chromosomes/contigs shared between the fasta and GTF file you have selected!\n")


# Default memory allocation: 10G
rule bowtie2Index:
    input: genome_fasta
    output: os.path.join(outdir, "BowtieIndex/genome.rev.2.bt2")
    params:
      basedir = os.path.join(outdir, "BowtieIndex")
    conda: CONDA_CREATE_INDEX_ENV
    threads: lambda wildcards: 10 if 10<max_thread else max_thread
    shell: """
        ln -s {input} {params.basedir}/genome.fa
        bowtie2-build -t {threads} {params.basedir}/genome.fa {params.basedir}/genome
        if [[ -f BowtieIndex/genome.rev.2.bt2l ]]; then ln -s genome.rev.2.bt2l {output} ; fi
        """

# Default memory allocation: 20G
rule hisat2Index:
    input: genome_fasta
    output: os.path.join(outdir, "HISAT2Index/genome.6.ht2")
    params:
      basedir = os.path.join(outdir, "HISAT2Index")
    threads: lambda wildcards: 10 if 10<max_thread else max_thread
    conda: CONDA_CREATE_INDEX_ENV
    shell: """
        ln -s {input} {params.basedir}/genome.fa
        hisat2-build -q -p {threads} {params.basedir}/genome.fa {params.basedir}/genome
        """


# Default memory allocation: 1G
rule makeKnownSpliceSites:
    input: genes_gtf
    output: known_splicesites
    conda: CONDA_CREATE_INDEX_ENV
    threads: lambda wildcards: 10 if 10<max_thread else max_thread
    shell: """
        hisat2_extract_splice_sites.py {input} > {output}
        """


# Default memory allocation: 80G
rule starIndex:
    input: genome_fasta
    output: os.path.join(outdir, "STARIndex/SAindex")
    params:
      basedir = os.path.join(outdir, "STARIndex")
    conda: CONDA_CREATE_INDEX_ENV
    threads: lambda wildcards: 10 if 10<max_thread else max_thread
    shell: """
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.basedir} --genomeFastaFiles {input}
        if [[ -w Log.out ]]; then rm -v Log.out; elif [[ -w {params.basedir}/Log.out ]]; then rm -v {params.basedir}/Log.out; fi
        """

rule genes_bed2fasta:
    input:
        bed = genes_bed,
        genome_fasta = genome_fasta
    output:
        "annotation/genes.fa"
    benchmark:
        "annotation/.benchmark/annotation_bed2fasta.benchmark"
    threads: 1
    conda: CONDA_CREATE_INDEX_ENV
    shell:
        "bedtools getfasta -name -s -split -fi {input.genome_fasta} -bed <(cat {input.bed} | cut -f1-12) | sed 's/(.*)//g' | sed 's/:.*//g' > {output}"


rule salmonIndex:
    input:
        "annotation/genes.fa",
        genome_fasta
    output:
        os.path.join(outdir, "SalmonIndex/decoys.txt"),
        temp(os.path.join(outdir, "SalmonIndex/seq.fa")),
        os.path.join(outdir, "SalmonIndex/seq.bin")
    params:
        salmonIndexOptions = salmonIndexOptions if salmonIndexOptions else ""
    threads: lambda wildcards: 16 if 16<max_thread else max_thread
    conda: CONDA_CREATE_INDEX_ENV
    shell: """
        grep "^>" {input[1]} | cut -d " " -f 1 | tr -d ">" > {output[0]}
        cat {input[0]} {input[1]} > {output[1]}
        salmon index -p {threads} -t {output[1]} -d {output[0]} -i SalmonIndex {params.salmonIndexOptions}
        """


##### the code for obtaining spliced/unspliced counts from Alevin is based on Soneson et al.2020, bioRxiv, https://doi.org/10.1101/2020.03.13.990069

rule run_eisaR:
    input:
        gtf = genes_gtf,
        genome_fasta = genome_fasta
    output:
        joint_fasta = temp(os.path.join(outdir, "annotation/cDNA_introns.joint.fa")),
        joint_t2g = os.path.join(outdir, "annotation/cDNA_introns.joint.t2g")
    params:
        wdir = os.path.join(outdir, "annotation"),
        scriptdir = workflow_rscripts,
        isoform_action = "separate",
        flank_length = eisaR_flank_length,
        gtf = lambda wildcards,input: os.path.join(outdir, input.gtf),
        joint_fasta = lambda wildcards,output: output.joint_fasta,
        joint_t2g = lambda wildcards,output: output.joint_t2g
    conda: CONDA_eisaR_ENV
    script: "../rscripts/scRNAseq_eisaR.R"



#uses decoys generated by rule SalmonIndex in Salmon.snakefile

rule Salmon_index_joint_fa:
    input:
        joint_fasta = os.path.join(outdir, "annotation/cDNA_introns.joint.fa"),
        decoys = os.path.join(salmon_index, "decoys.txt"),
        genome_fasta = genome_fasta
    output:
        seq_fa = temp(os.path.join(outdir, "SalmonIndex_RNAVelocity/seq.fa")),
        velo_index = os.path.join(outdir, "SalmonIndex_RNAVelocity/seq.bin")
    params:
        salmonIndexOptions = salmonIndexOptions
    threads: lambda wildcards: 16 if 16<max_thread else max_thread
    conda: CONDA_SALMON_ENV
    shell:"""
        cat {input.joint_fasta} {input.genome_fasta} > {output.seq_fa}
        salmon index -p {threads} -t {output.seq_fa} -d {input.decoys} -i SalmonIndex_RNAVelocity {params.salmonIndexOptions}
        """



# Default memory allocation: 8G
rule bwaIndex:
    input: genome_fasta
    output: os.path.join(outdir, "BWAIndex/genome.fa.sa")
    params:
      genome = os.path.join(outdir, "BWAIndex", "genome.fa")
    conda: CONDA_CREATE_INDEX_ENV
    shell: """
        ln -s {input} {params.genome}
        bwa index {params.genome}
        """

# Default memory allocation: 8G
rule bwamem2Index:
    input: genome_fasta
    output: os.path.join(outdir, "BWA-MEM2Index/genome.fa.bwt.2bit.64")
    params:
      genome = os.path.join(outdir, "BWA-MEM2Index", "genome.fa")
    conda: CONDA_CREATE_INDEX_ENV
    shell: """
        ln -s {input} {params.genome}
        bwa-mem2 index {params.genome}
        """


# Default memory allocation: 8G
rule bwamethIndex:
    input: genome_fasta
    output: os.path.join(outdir, "BWAmethIndex/genome.fa.bwameth.c2t.sa")
    params:
      genome = os.path.join(outdir, "BWAmethIndex", "genome.fa")
    conda: CONDA_CREATE_INDEX_ENV
    shell: """
        ln -s {input[0]} {params.genome}
        bwameth.py index {params.genome}
        """

# Default memory allocation: 8G
rule bwameth2Index:
    input: genome_fasta
    output: os.path.join(outdir, "BWAmeth2Index/genome.fa.bwameth.c2t.bwt.2bit.64")
    params:
      genome = os.path.join(outdir, "BWAmeth2Index", "genome.fa")
    conda: CONDA_CREATE_INDEX_ENV
    shell: """
        ln -s {input[0]} {params.genome}
        bwameth.py index-mem2 {params.genome}
        """

# Default memory allocation: 1G
rule copyBlacklist:
    output: os.path.join(outdir, "annotation/blacklist.bed")
    params:
        url = blacklist
    run:
        downloadFile(params.url, output)

rule copySpikeinBlacklist:
    output: temp(os.path.join(outdir, "annotation/spikein.blacklist_ori.bed"))
    params:
        url = spikeinBlacklist
    run:
        downloadFile(params.url, output)

rule renameSpikeinChromsBlacklist:
    input:  os.path.join(outdir,"annotation/spikein.blacklist_ori.bed")
    output: spikein_blacklist_bed
    params:
        spikeinExt = spikeinExt
    shell: """
        awk -v FS='\\t' -v OFS='\\t' '{{ if($1 !~ /^#/){{$1=$1\"{params.spikeinExt}\"; print $0}} else{{print $0}} }}' {input} > {output}
    """


# Default memory allocation: 1G
rule computeEffectiveGenomeSize:
    input: genome_fasta if not spikeinGenomeURL else os.path.join(outdir,"genome_fasta/host.genome.fa")
    output: os.path.join(outdir, "genome_fasta", "effectiveSize")
    conda: CONDA_SHARED_ENV
    shell: """
        seqtk comp {input} | awk '{{tot += $3 + $4 + $5 + $6}}END{{print tot}}' > {output}
        """
