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
            sed 's/\s\+/{params.spikeinExt} /' {input} > {output}
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

if tesmall:
    # Ported from smallrna_pipeline_snakemake/Snakefile rather than shelling out to
    # its standalone build_smallrna_reference.py script, so it plugs into this DAG
    # (reuses genome_fasta/genes_gtf, participates in -n/--dag, etc) instead of being
    # an opaque black-box step. The stage-module functions it calls are vendored at
    # shared/tools/smallrna/ (imported up in workflows/createIndices/Snakefile).
    #
    # tesmallGenome (defaults to `genome`, overridable via --tesmallGenome) drives
    # RepeatMasker/miRBase/piRNAdb lookups. It's validated against common.py's
    # normalize_genome() up in createIndices.py before snakemake even starts, since
    # `genome` here is just a free-form label for the output YAML's filename, not
    # necessarily a real build code.
    #
    # rmsk.txt is fetched fresh here rather than reusing rmsk_file (from --rmskURL):
    # rmsk_file is whatever raw format the user's URL happens to serve, while
    # extract_te_from_rmsk/extract_structural_rna_from_rmsk require the specific
    # Ensembl-coordinate conversion fetch_rmsk_txt() does (UCSC download, 'chr'
    # stripped, 1-based start). Reusing rmsk_file directly would silently give
    # wrong coordinates whenever its format doesn't happen to match that.

    # TEsmall/genomes/<tesmallGenome>/{annotation,sequence} -- nested under the
    # build name so multiple genome versions can eventually coexist under one
    # TEsmall/ root without collisions (each build gets its own rmsk.txt,
    # .chains/, .build_info/, sequence/, annotation/).
    TESMALL_GENOME = tesmallGenome
    TESMALL_RMSK_ASSEMBLY = normalize_genome(TESMALL_GENOME)

    TESMALL_OUT = os.path.join(outdir, "TEsmall")
    TESMALL_GENOME_DIR = os.path.join(TESMALL_OUT, "genomes", TESMALL_GENOME)
    TESMALL_SEQ_DIR = os.path.join(TESMALL_GENOME_DIR, "sequence")
    TESMALL_ANNOTATION_DIR = os.path.join(TESMALL_GENOME_DIR, "annotation")
    TESMALL_CHAIN_DIR = os.path.join(TESMALL_GENOME_DIR, ".chains")
    TESMALL_BUILD_INFO_DIR = os.path.join(TESMALL_GENOME_DIR, ".build_info")

    TESMALL_RMSK_TXT = os.path.join(TESMALL_ANNOTATION_DIR, "rmsk.txt")
    TESMALL_EXON_BED = os.path.join(TESMALL_ANNOTATION_DIR, "exon.bed")
    TESMALL_INTRON_BED = os.path.join(TESMALL_ANNOTATION_DIR, "intron.bed")
    TESMALL_STRUCTURAL_BED = os.path.join(TESMALL_ANNOTATION_DIR, "structural_RNA.bed")
    TESMALL_TE_BED = os.path.join(TESMALL_ANNOTATION_DIR, "TE.bed")
    TESMALL_TDNA_FA = os.path.join(TESMALL_SEQ_DIR, "tDNA.fa")
    TESMALL_RDNA_FA = os.path.join(TESMALL_SEQ_DIR, "rDNA.fa")
    TESMALL_GENOME_FA = os.path.join(TESMALL_SEQ_DIR, "genome.fa")
    TESMALL_HAIRPIN_BED = os.path.join(TESMALL_ANNOTATION_DIR, "hairpin.bed")
    TESMALL_MATURE_BED = os.path.join(TESMALL_ANNOTATION_DIR, "miRNA.bed")
    TESMALL_PIRNA_BED = os.path.join(TESMALL_ANNOTATION_DIR, "piRNA_cluster.bed")
    TESMALL_BOWTIE_DIR = os.path.join(TESMALL_SEQ_DIR, "bowtie_index")
    TESMALL_PREFIX_FOR_NAME = {"tDNA": "sncRNA:tRNA:", "rDNA": "sncRNA:rRNA:"}
    # TEsmall (the analysis tool) expects bowtie1 indexes of genome/tDNA/rDNA inside
    # its --dbfolder, so unlike smallrna_pipeline_snakemake's own build_bowtie
    # config toggle, these are always built here -- there's no point in a
    # --tesmall run that doesn't produce something TEsmall can actually load.
    TESMALL_FASTA_FOR_BOWTIE = {"genome": TESMALL_GENOME_FA, "tDNA": TESMALL_TDNA_FA, "rDNA": TESMALL_RDNA_FA}

    rule tesmall_link_genome:
        input:
            fa=genome_fasta,
            fai=genome_index,
        output:
            fa=TESMALL_GENOME_FA,
            fai=TESMALL_GENOME_FA + ".fai",
        run:
            link_genome_fasta(os.path.join(outdir, "genome_fasta"), TESMALL_SEQ_DIR)

    rule tesmall_fetch_rmsk:
        output:
            TESMALL_RMSK_TXT,
        run:
            ok = fetch_rmsk_txt(TESMALL_RMSK_ASSEMBLY, output[0])
            if not ok:
                raise RuntimeError(f"Failed to download RepeatMasker table for {TESMALL_RMSK_ASSEMBLY}")

    rule tesmall_exon_intron_raw:
        input:
            gtf=genes_gtf,
        output:
            exon=os.path.join(TESMALL_ANNOTATION_DIR, "exon.raw.bed"),
            intron=os.path.join(TESMALL_ANNOTATION_DIR, "intron.raw.bed"),
        run:
            write_raw_exon_bed(input.gtf, output.exon)
            write_raw_intron_bed(input.gtf, output.intron)

    rule tesmall_structural_rna_raw:
        input:
            gtf=genes_gtf,
            rmsk=TESMALL_RMSK_TXT,
        output:
            os.path.join(TESMALL_ANNOTATION_DIR, "structural_RNA.raw.bed"),
        run:
            write_raw_structural_rna(input.gtf, input.rmsk, output[0])

    rule tesmall_te_bed:
        input:
            rmsk=TESMALL_RMSK_TXT,
        output:
            TESMALL_TE_BED,
        run:
            extract_te_from_rmsk(input.rmsk, output[0])

    rule tesmall_collapse_bed:
        input:
            raw=os.path.join(TESMALL_ANNOTATION_DIR, "{name}.raw.bed"),
        output:
            os.path.join(TESMALL_ANNOTATION_DIR, "{name}.bed"),
        wildcard_constraints:
            name="exon|intron|structural_RNA",
        conda: CONDA_CREATE_TEsmall_ENV
        shell:
            r"""
            bedtools sort -i {input.raw} | \
            bedtools groupby -g 1,2,3,6 -c 4 -o collapse | \
            awk -F'\t' -v OFS='\t' '{{print $1,$2,$3,$5,0,$4}}' > {output}
            rm -f {input.raw}
            """

    rule tesmall_select_structural_subset:
        input:
            bed=TESMALL_STRUCTURAL_BED,
        output:
            os.path.join(TESMALL_SEQ_DIR, "{name}.selected.bed"),
        wildcard_constraints:
            name="tDNA|rDNA",
        params:
            prefix=lambda wc: TESMALL_PREFIX_FOR_NAME[wc.name],
        run:
            filter_bed_by_prefix(input.bed, params.prefix, output[0])

    rule tesmall_extract_fasta:
        input:
            genome_fa=TESMALL_GENOME_FA,
            bed=os.path.join(TESMALL_SEQ_DIR, "{name}.selected.bed"),
        output:
            os.path.join(TESMALL_SEQ_DIR, "{name}.raw.fa"),
        wildcard_constraints:
            name="tDNA|rDNA",
        conda: CONDA_CREATE_TEsmall_ENV
        shell:
            r"""
            if [ -s {input.bed} ]; then
                bedtools getfasta -s -name -fi {input.genome_fa} -bed {input.bed} -fo {output}
            else
                : > {output}
            fi
            rm -f {input.bed}
            """

    rule tesmall_finalize_fasta:
        input:
            raw=os.path.join(TESMALL_SEQ_DIR, "{name}.raw.fa"),
        output:
            os.path.join(TESMALL_SEQ_DIR, "{name}.fa"),
        wildcard_constraints:
            name="tDNA|rDNA",
        run:
            if wildcards.name == "tDNA":
                fix_trna_headers(input.raw, output[0])
            else:
                shutil.copyfile(input.raw, output[0])
            os.remove(input.raw)

    rule tesmall_faidx:
        input:
            fa=os.path.join(TESMALL_SEQ_DIR, "{name}.fa"),
        output:
            os.path.join(TESMALL_SEQ_DIR, "{name}.fa.fai"),
        wildcard_constraints:
            name="tDNA|rDNA",
        conda: CONDA_CREATE_TEsmall_ENV
        shell:
            r"""
            if [ -s {input.fa} ]; then
                samtools faidx {input.fa}
            else
                : > {output}
            fi
            """

    checkpoint tesmall_mirna_raw:
        output:
            hairpin=os.path.join(TESMALL_ANNOTATION_DIR, "hairpin.raw.bed"),
            mature=os.path.join(TESMALL_ANNOTATION_DIR, "miRNA.raw.bed"),
            build_info=os.path.join(TESMALL_BUILD_INFO_DIR, "mirna_build_info.txt"),
        run:
            ensure_dir(TESMALL_BUILD_INFO_DIR)
            download_smallrna_raw("mirbase", TESMALL_GENOME, TESMALL_ANNOTATION_DIR, build_info_name="mirna_build_info.txt")
            shutil.move(os.path.join(TESMALL_ANNOTATION_DIR, "mirna_build_info.txt"), output.build_info)

    checkpoint tesmall_pirna_raw:
        output:
            pirna=os.path.join(TESMALL_ANNOTATION_DIR, "piRNA_cluster.raw.bed"),
            build_info=os.path.join(TESMALL_BUILD_INFO_DIR, "pirna_build_info.txt"),
        run:
            ensure_dir(TESMALL_BUILD_INFO_DIR)
            download_smallrna_raw("pirnadb", TESMALL_GENOME, TESMALL_ANNOTATION_DIR, build_info_name="pirna_build_info.txt")
            shutil.move(os.path.join(TESMALL_ANNOTATION_DIR, "pirna_build_info.txt"), output.build_info)

    def _tesmall_read_build_info(path):
        src, tgt = open(path).read().split()
        return src, tgt

    def tesmall_mirna_finalize_input(wildcards):
        build_info = checkpoints.tesmall_mirna_raw.get().output.build_info
        src, tgt = _tesmall_read_build_info(build_info)
        suffix = "raw" if src == tgt else "lifted"
        return {
            "hairpin": os.path.join(TESMALL_ANNOTATION_DIR, f"hairpin.{suffix}.bed"),
            "mature": os.path.join(TESMALL_ANNOTATION_DIR, f"miRNA.{suffix}.bed"),
        }

    def tesmall_pirna_finalize_input(wildcards):
        build_info = checkpoints.tesmall_pirna_raw.get().output.build_info
        src, tgt = _tesmall_read_build_info(build_info)
        suffix = "raw" if src == tgt else "lifted"
        return {"pirna": os.path.join(TESMALL_ANNOTATION_DIR, f"piRNA_cluster.{suffix}.bed")}

    rule tesmall_mirna_finalize:
        input:
            unpack(tesmall_mirna_finalize_input),
        output:
            hairpin=TESMALL_HAIRPIN_BED,
            mature=TESMALL_MATURE_BED,
        run:
            shutil.copyfile(input.hairpin, output.hairpin)
            shutil.copyfile(input.mature, output.mature)
            for f in (input.hairpin, input.mature):
                if os.path.exists(f):
                    os.remove(f)
            raw = checkpoints.tesmall_mirna_raw.get().output
            for f in (raw.hairpin, raw.mature):
                if os.path.exists(f):
                    os.remove(f)

    rule tesmall_pirna_finalize:
        input:
            unpack(tesmall_pirna_finalize_input),
        output:
            TESMALL_PIRNA_BED,
        run:
            shutil.copyfile(input.pirna, output[0])
            if os.path.exists(input.pirna):
                os.remove(input.pirna)
            raw = checkpoints.tesmall_pirna_raw.get().output.pirna
            if os.path.exists(raw):
                os.remove(raw)

    _TESMALL_MIRNA_NAMES = {"hairpin", "miRNA"}

    def tesmall_chain_for(wildcards):
        if wildcards.name in _TESMALL_MIRNA_NAMES:
            build_info = checkpoints.tesmall_mirna_raw.get().output.build_info
        else:
            build_info = checkpoints.tesmall_pirna_raw.get().output.build_info
        src, tgt = _tesmall_read_build_info(build_info)
        chain_name = f"{src}To{tgt[0].upper()}{tgt[1:]}.over.chain.gz"
        return os.path.join(TESMALL_CHAIN_DIR, chain_name)

    rule tesmall_download_chain:
        output:
            os.path.join(TESMALL_CHAIN_DIR, "{src}To{tgt_cap}.over.chain.gz"),
        run:
            ensure_dir(TESMALL_CHAIN_DIR)
            url = f"{UCSC_BASE}/{wildcards.src}/liftOver/{wildcards.src}To{wildcards.tgt_cap}.over.chain.gz"
            download_or_raise(url, output[0])

    rule tesmall_liftover_bed:
        input:
            raw=os.path.join(TESMALL_ANNOTATION_DIR, "{name}.raw.bed"),
            chain=tesmall_chain_for,
        output:
            os.path.join(TESMALL_ANNOTATION_DIR, "{name}.lifted.bed"),
        wildcard_constraints:
            name="hairpin|miRNA|piRNA_cluster",
        conda: CONDA_CREATE_TEsmall_ENV
        shell:
            r"""
            awk -F'\t' -v OFS='\t' '{{
                c=$1
                if (c !~ /^chr/) {{ if (c=="MT" || c=="mt") c="chrM"; else c="chr" c }}
                $1=c; print
            }}' {input.raw} > {input.raw}.chr.bed

            liftOver {input.raw}.chr.bed {input.chain} {output}.chr.lifted {output}.unmapped

            n=$(grep -vc '^#' {output}.unmapped 2>/dev/null || true)
            if [ -n "$n" ] && [ "$n" -gt 0 ]; then
                echo "[WARN] liftOver {wildcards.name}: $n interval(s) dropped" >&2
            fi

            awk -F'\t' -v OFS='\t' '{{
                c=$1
                if (c ~ /^chr/) {{ c=substr(c,4); if (c=="M") c="MT" }}
                $1=c; print
            }}' {output}.chr.lifted > {output}

            rm -f {input.raw}.chr.bed {output}.chr.lifted {output}.unmapped
            """

    rule tesmall_bowtie_index:
        input:
            fasta=lambda wc: TESMALL_FASTA_FOR_BOWTIE[wc.name],
        output:
            multiext(os.path.join(TESMALL_BOWTIE_DIR, "{name}"),
                      ".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt",
                      ".rev.1.ebwt", ".rev.2.ebwt"),
        wildcard_constraints:
            name="genome|tDNA|rDNA",
        params:
            prefix=lambda wc: os.path.join(TESMALL_BOWTIE_DIR, wc.name),
        threads: lambda wildcards: 4 if 4 < max_thread else max_thread
        conda: CONDA_CREATE_TEsmall_ENV
        shell:
            "bowtie-build --threads {threads} {input.fasta} {params.prefix}"

    rule create_tesmall:
        input:
            TESMALL_EXON_BED, TESMALL_INTRON_BED, TESMALL_STRUCTURAL_BED, TESMALL_TE_BED,
            TESMALL_HAIRPIN_BED, TESMALL_MATURE_BED, TESMALL_PIRNA_BED,
            os.path.join(TESMALL_SEQ_DIR, "tDNA.fa.fai"), os.path.join(TESMALL_SEQ_DIR, "rDNA.fa.fai"),
            TESMALL_GENOME_FA + ".fai",
            [multiext(os.path.join(TESMALL_BOWTIE_DIR, name),
                      ".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt")
             for name in TESMALL_FASTA_FOR_BOWTIE],
        output: touch(os.path.join(TESMALL_GENOME_DIR, ".tesmall_done"))
        run:
            # rmsk.txt is a pure intermediate for structural_RNA.bed/TE.bed above;
            # remove it now that both have consumed it. .chains/ is likewise only
            # needed transiently by liftover_bed. .build_info/ is NOT touched here
            # -- deleting a checkpoint's own output risks Snakemake silently
            # rerunning the whole checkpoint (re-downloading from miRBase/piRNAdb)
            # the next time something re-resolves it (snakemake#609), so it just
            # stays in its own hidden dir instead of ever needing to be cleaned up.
            if os.path.exists(TESMALL_RMSK_TXT):
                os.remove(TESMALL_RMSK_TXT)
            if os.path.isdir(TESMALL_CHAIN_DIR):
                shutil.rmtree(TESMALL_CHAIN_DIR)
