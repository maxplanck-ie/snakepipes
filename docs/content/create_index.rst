Creating indices
================

While snakePipes makes it convenient to install packages, it does not yet make creating the genome and other index files these packages sometimes need convenient (this is planned but not yet implemented). The location and structure of the indices provided in the organism YAML files supplied with snakePipes matches the `structure we use internally <https://github.com/maxplanck-ie/data_repository/tree/master/organisms>`__. As modifying Makefiles is not a task anyone looks forward too, we will list below the basics of what is needed to generate the various indices used by snakePipes. For convenience, each section will be titled with the variable or variables in the YAML file that it refers to. We will be using the human genome as an example here.

.. note:: The ordering below is intentional, as some files are needed to create others.


genome_fasta
------------

Everything starts with a full-genome fasta file. We typically download these from Ensembl, which provides a very convenient `download page <http://www.ensembl.org/info/data/ftp/index.html>`__. Care should be made when you decide what fasta file to download. Most importantly, the source of your fasta file (e.g., Ensembl) should also provide annotations in GTF format. Mixing and matching fasta and GTF files from different sources is a recipe for disaster! Further, always make sure to download the DNA (not cDNA, CDS, etc.) fasta file. Some providers, such as Ensembl, provide "sm" (soft-masked) and "rm" (hard-masked) files as well as "primary" and "toplevel" files. In such cases, download the "sm toplevel" fasta file.

The exact commands to then download and prepare the file are::

    wget ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
    gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
    mv Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa genome.fa

You can then put the location of `genome.fa` in the `genome_fasta` field.


genome_index
------------

This is created from the file listed in `genome_fasta`. Assuming that is set to `genome.fa` you can produce the requisite file with::

    conda install -c bioconda -c conda-forge samtools
    samtools index genome.fa

This will then produce a file ending in ".fai" that you can put in the `genome_index` field.


genome_2bit
-----------

The 2bit format is a randomly accessible representation of a genome. Within snakePipes, this is only used if you specify the `--gcbias` option in the DNA-mapping pipeline. To create a 2bit file, you first need the file specified in `genome_fasta`::

    conda install -c bioconda -c conda-forge ucsc-fatotwobit
    faToTwoBit genome.fa genome.2bit

You can then put the location of `genome.2bit` in the `genome_2bit` field.


genome_size
-----------

A number of commands need to know the "effective genome size", such as `bamCoverage` and `MACS2`. For heavily used genomes you can simply search google for a reasonable value. For less commonly used genomes, you can approximate this be counting the total number of bases in the genome and subtracting the number of Ns. There are MANY ways to do this, some of which are `listed here <https://www.biostars.org/p/19426/>`__. When in doubt, simply put the genome size here (you'll find it on the last line of the file listed in `genome_index`).


genes_gtf
---------

Many tools in snakePipes need to know where genes are located and their structure. You can typically download GTF files from Ensembl or other sources. For the human genome the following will suffice::

    wget ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.chr.gtf.gz
    gunzip Homo_sapiens.GRCh38.93.chr.gtf.gz
    mv Homo_sapiens.GRCh38.93.chr.gtf genes.gtf

You can then put the location of `genes.gtf` in the `genes_gtf` field.


genes_bed
---------

The BED format is an alternative representation of much of the information in a GTF file. This is needed by a number of tools in snakePipes, since it's easier for tools to parse than a GTF file. Given the file listed in `genes_gtf`, you can create this file as follows::

    conda install -c bioconda -c conda-forge ucsc-gtftogenepred ucsc-genepredtobed
    awk '{if ($$3 != "gene") print $$0;}' genes.gtf \
        | grep -v "^#" \
        | gtfToGenePred /dev/stdin /dev/stdout \
        | genePredToBed stdin genes.bed

You can then place the location of `genes.bed` in the `genes_bed` field.


bowtie2_index
-------------

Bowtie2 is used in the DNA-mapping pipeline. It requires a custom index of the file listed in `genome_fasta`. It's often most convenient to create a new directory for this index and to create a symlink in it to the file in `genome_fasta`. Assuming the file in `genome_fasta` is located in `genome_fasta/genome.fa`::

    conda install -c bioconda -c conda-forge bowtie2
    mkdir BowtieIndex
    ln -s genome_fasta/genome.fa BowtieIndex/genome.fa
    bowtie2-build BowtieIndex/genome.fa BowtieIndex/genome

You can then place `/path/to/BowtieIndex/genome` in the `bowtie2_index` field. Note that there will not be an actual file with this name, rather it serves as a prefix for all of the files that Bowtie2 needs.


hisat2_index
------------

The RNA-seq pipeline can be instructed to use HISAT2 for alignment rather than STAR. As with Bowtie2, this tool requires its own index, which can be most conveniently placed in a new directory and is made from the file listed in `genome_fasta`::

    conda install -c bioconda -c conda-forge hisat2
    mkdir HISAT2Index
    ln -s genome_fasta/genome.fa HISAT2Index/genome.fa
    hisat2-build -q -p 10 HISAT2Index/genome.fa BowtieIndex/genome

This uses 10 threads for indexing, please change that as appropriate. You can then place `/path/to/HISAT2Index/genome` in the `hisat2_index` field. Note that there will not be an actual file with this name, rather it serves as a prefix for all of the files that HISAT2 needs.


known_splicesites
-----------------

The HISAT2 index can be created with a particular genome annotation. In our experience, it's most convenient to simply index the genome and then later add the splice sites from a given annotation into an index at run time. This facilitates updating versions of GTF files without needing to reindex the genome each time. Thus, snakePipes expects a file of splice sites to give to HISAT2. This requires the file listed in `genes_gtf`::

    conda install -c bioconda -c conda-forge hisat2
    mkdir HISAT2Index
    extract_splice_sites.py genes.gtf > HISAT2Index/splice_sites.txt

You can then place the location of `HISAT2Index/splice_sites.txt` in the `known_splicesites` field.


star_index
----------

The RNA-seq pipeline uses STAR for alignment by default. As with other aligners, this requires a custom index to be made from the file listed in `genome_fasta`::

    conda install -c bioconda -c conda-forge star
    mkdir STARIndex
    STAR --runThreadN 10 --runMode genomeGenerate --genomeDir STARIndex --genomeFastaFiles genome_fasta/genome.fa
    rm -f Log.out

You can then place the location of the `STARIndex` directory in the `star_index` field.


bwa_index
---------

BWA is used for alignment in the Hi-C pipeline. To create an index for BWA, you need the file listed in `genome_fasta`::

    conda install -c bioconda -c conda-forge bwa
    mkdir BWAindex
    ln -s genome_fasta/genome.fa BWAindex/genome.fa
    bwa index BWAindex/genome.fa

You can then place the location of `BWAindex/genome.fa` in the `bwa_index` field.


bwameth_index
-------------

The WGBS pipeline uses BWA-meth for alignment. To create the index for it, you need the file listed in `genome_fasta`::


    conda install -c bioconda -c conda-forge bwameth
    mkdir BWAmethIndex
    ln -s genome_fasta/genome.fa BWAmethIndex/genome.fa
    bwameth.py index BWAmethIndex/genome.fa

You can then place the location of `BWAmethIndex/genome.fa` in the `bwameth_index` field.


extended_coding_regions_gtf
---------------------------

The RNA-seq pipeline can be instructed to compute the level of genomic contamination, as done by some consortia. This requires the file from `genes_gtf`::

    grep -v "^#" gene.gtf \
        | awk 'BEGIN{OFS="\t"}{if($$3 == "gene" || $$3 == "transcript") {$$4 -= 500; $$5 += 500; if($$4 < 1) {$$4 = 1}; print}}' > genes.slop.gtf

You can then place the location of `genes.slop.gtf` in the `extended_coding_regions_gtf` field.


blacklist_bed
-------------

Some protocols, in particular ChIP-seq, produce spurious signal in the same regions regardless of experiment. These regions then need to be "blacklisted" or ignored in the analysis. Such files are typically produced by consortia such as ENCODE::

    wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz
    gunzip hg38.blacklist.bed.gz

You can then place the location of `hg38.blacklist.bed` in the `blacklist_bed` field.

.. note:: Blacklists only exist for very heavily used genomes, such as mouse and human. Feel free to leave this blank if no such file exists for your genome.

ignore_forNorm
--------------

A number of chromosomes can skew normalization. The classic example of this is sex chromsomes and mitochondria, which will vary between individuals. Many programs in DeepTools can be instructed to ignore these by placing them as a space-separated list in the `ignore_forNorm` field. We additionally find it convenient to exclude unplaced contigs.
