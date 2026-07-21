Cookbook
========

.. contents::
   :local:
   :depth: 2

--------------------------------------------------------------------
Create Hybrid Genome with Spike-in Sequences (:ref:`createIndices`)
--------------------------------------------------------------------

Generate a reference genome including both the native organism (referred to as host)  and spike-in control sequences.   

.. code-block:: bash

   createIndices -o customIndices/GRch38_dm6 --tools bowtie2  --genomeURL organisms/GRCh38_ensembl/genome_fasta/genome.fa --gtfURL organisms/GRCh38_ensembl/gencode/release_31/genes.gtf --blacklist organisms/GRCh38_ensembl/akundaje/blacklist.UseMe.bed --spikeinGenomeURL organisms/dm6_ensembl/genome_fasta/genome.fa --spikeinGtfURL organisms/dm6_ensembl/ensembl/release-96/genes.gtf customIndices/GRch38_dm6/GRCh38_g31_dm6


-------------------------------------------------------------
Cut&Tag Data Analysis (:ref:`DNAmapping` and :ref:`ChIPseq`)
-------------------------------------------------------------

Process and analyze Cut&Tag data for chromatin binding. Fastq files are mapped to a host-spikein hybrid genome and spikein sequences are used for normalization of bam coverage tracks.   

.. code-block:: bash

   DNAmapping --cutntag --trim --trimmerOptions ' -a nexteraF=CTGTCTCTTATA -A nexteraR=CTGTCTCTTATA ' --fastqc --dedup --mapq 3 -i $input_folder -o analysis_dedup customIndices/GRch38_dm6/GRCh38_g31_dm6.yaml

.. code-block:: bash

   ChIPseq -d analysis_ChIPseq --fromBAM analysis_dedup/filtered_bam --bamExt .filtered.bam --cutntag customIndices/GRch38_dm6/GRCh38_g31_dm6.yaml chip_seq_sample_config.yaml


--------------------------------------------------------
Differential Binding on Target Regions (:ref:`ChIPseq`)
--------------------------------------------------------

Identify regions showing significant changes in  binding between conditions.

.. code-block:: bash

   ChIPseq -d analysis_chipseq --fromBAM analysis_dna --bamExt .bam  --externalBed rmsk.bed  --sampleSheet sampleSheet.tsv mm10_gencodeM19 chip_seq_sample_config.yaml



------------------------------------------------------
Time-course Analysis of mRNAseq Data (:ref:`mRNAseq`)
------------------------------------------------------

Investigate temporal changes in transcript levels across multiple time points by leveraging an LRT test. For this implementation, the `condition` column in the sample sheet should contain groups corresponding to time points.

.. code-block:: bash

   mRNAseq -i RNAseq -o analysis --sampleSheet sampleSheet.csv  --LRT mm10_gencodeM19 


-----------------------------------------------------------------------------
Differential Expression of Transcribed Repetitive Elements (:ref:`ncRNAseq`)
-----------------------------------------------------------------------------

Assess changes in expression for transcribed repeats. The organism yaml must contain the key rmsk_file pointing to the repeat masker txt file.

.. code-block:: bash

   ncRNAseq -i RNAseq -o analysis --sampleSheet sampleSheet.csv  mm10_gencodeM19


---------------------------------------------------
Hi-C Data Analysis with TAD Calling (:ref:`HiC`)
---------------------------------------------------

Analyze chromatin conformation and identify topologically associating domains (TADs).

.. code-block:: bash

   HiC -i merged_fq -o analysis --fastqc --trim --enzyme DpnII --binSize 5000 mm10_gencodeM19


------------------------------------------------------
Allele-specific Hi-C Data Analysis (:ref:`makePairs`)
------------------------------------------------------

Distinguish allelic differences in chromatin interactions.

.. code-block:: bash

   makePairs -i input-dir -o output-dir --VCFfile vcf --strains s1,s2 dm6


-------------------------------------------------------------
Differential Transcript Expression Analysis (:ref:`mRNAseq`)
-------------------------------------------------------------

Quantify and compare transcript expression between different conditions.

.. code-block:: bash

   mRNAseq -i RNAseq -o analysis -m alignment-free --sampleSheet sampleSheet.csv mm10_gencodeM19

