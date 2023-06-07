#!/bin/bash
set -ex
if [[ ${CI:-"false"} == "true" ]]; then
    export PATH="$HOME/miniconda/bin:$PATH"
    hash -r
    python -m pip install --no-deps --ignore-installed .
fi

# Needed by DNA, HiC, mRNA-seq, WGBS and scRNA-seq workflows
mkdir -p PE_input
touch PE_input/sample1_R1.fastq.gz PE_input/sample1_R2.fastq.gz \
      PE_input/sample2_R1.fastq.gz PE_input/sample2_R2.fastq.gz \
      PE_input/sample3_R1.fastq.gz PE_input/sample3_R2.fastq.gz \
      PE_input/sample4_R1.fastq.gz PE_input/sample4_R2.fastq.gz \
      PE_input/sample5_R1.fastq.gz PE_input/sample5_R2.fastq.gz \
      PE_input/sample6_R1.fastq.gz PE_input/sample6_R2.fastq.gz
mkdir -p SE_input
touch SE_input/sample1_R1.fastq.gz \
      SE_input/sample2_R1.fastq.gz \
      SE_input/sample3_R1.fastq.gz \
      SE_input/sample4_R1.fastq.gz \
      SE_input/sample5_R1.fastq.gz \
      SE_input/sample6_R1.fastq.gz
# Needed by ChIP and ATAC workflows
mkdir -p BAM_input/deepTools_qc/bamPEFragmentSize BAM_input/filtered_bam BAM_input/Sambamba BAM_input/bamCoverage
touch BAM_input/sample1.bam \
      BAM_input/sample2.bam \
      BAM_input/sample3.bam \
      BAM_input/sample4.bam \
      BAM_input/sample5.bam \
      BAM_input/sample6.bam \
      BAM_input/filtered_bam/sample1.filtered.bam \
      BAM_input/filtered_bam/sample2.filtered.bam \
      BAM_input/filtered_bam/sample3.filtered.bam \
      BAM_input/filtered_bam/sample4.filtered.bam \
      BAM_input/filtered_bam/sample5.filtered.bam \
      BAM_input/filtered_bam/sample6.filtered.bam \
      BAM_input/filtered_bam/sample1.filtered.bam.bai \
      BAM_input/filtered_bam/sample2.filtered.bam.bai \
      BAM_input/filtered_bam/sample3.filtered.bam.bai \
      BAM_input/filtered_bam/sample4.filtered.bam.bai \
      BAM_input/filtered_bam/sample5.filtered.bam.bai \
      BAM_input/filtered_bam/sample6.filtered.bam.bai \
      BAM_input/Sambamba/sample1.markdup.txt \
      BAM_input/Sambamba/sample2.markdup.txt \
      BAM_input/Sambamba/sample3.markdup.txt \
      BAM_input/Sambamba/sample4.markdup.txt \
      BAM_input/Sambamba/sample5.markdup.txt \
      BAM_input/Sambamba/sample6.markdup.txt \
      BAM_input/deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv \
      BAM_input/bamCoverage/sample1.filtered.seq_depth_norm.bw \
      BAM_input/bamCoverage/sample2.filtered.seq_depth_norm.bw \
      BAM_input/bamCoverage/sample3.filtered.seq_depth_norm.bw \
      BAM_input/bamCoverage/sample4.filtered.seq_depth_norm.bw \
      BAM_input/bamCoverage/sample5.filtered.seq_depth_norm.bw \
      BAM_input/bamCoverage/sample6.filtered.seq_depth_norm.bw

mkdir -p allelic_BAM_input/allelic_bams allelic_BAM_input/filtered_bam  allelic_BAM_input/deepTools_qc/bamPEFragmentSize allelic_BAM_input/Sambamba allelic_BAM_input/bamCoverage/allele_specific
touch allelic_BAM_input/allelic_bams/sample1.genome1.sorted.bam \
      allelic_BAM_input/allelic_bams/sample1.genome2.sorted.bam \
      allelic_BAM_input/allelic_bams/sample2.genome1.sorted.bam \
      allelic_BAM_input/allelic_bams/sample2.genome2.sorted.bam \
      allelic_BAM_input/allelic_bams/sample3.genome1.sorted.bam \
      allelic_BAM_input/allelic_bams/sample3.genome2.sorted.bam \
      allelic_BAM_input/allelic_bams/sample4.genome1.sorted.bam \
      allelic_BAM_input/allelic_bams/sample4.genome2.sorted.bam \
      allelic_BAM_input/allelic_bams/sample5.genome1.sorted.bam \
      allelic_BAM_input/allelic_bams/sample5.genome2.sorted.bam \
      allelic_BAM_input/allelic_bams/sample6.genome1.sorted.bam \
      allelic_BAM_input/allelic_bams/sample6.genome2.sorted.bam \
      allelic_BAM_input/allelic_bams/sample1.genome1.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample1.genome2.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample2.genome1.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample2.genome2.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample3.genome1.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample3.genome2.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample4.genome1.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample4.genome2.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample5.genome1.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample5.genome2.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample6.genome1.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample6.genome2.sorted.bam.bai \
      allelic_BAM_input/deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv \
      allelic_BAM_input/filtered_bam/sample1.filtered.bam \
      allelic_BAM_input/filtered_bam/sample2.filtered.bam \
      allelic_BAM_input/filtered_bam/sample3.filtered.bam \
      allelic_BAM_input/filtered_bam/sample4.filtered.bam \
      allelic_BAM_input/filtered_bam/sample5.filtered.bam \
      allelic_BAM_input/filtered_bam/sample6.filtered.bam \
      allelic_BAM_input/filtered_bam/sample1.filtered.bam.bai \
      allelic_BAM_input/filtered_bam/sample2.filtered.bam.bai \
      allelic_BAM_input/filtered_bam/sample3.filtered.bam.bai \
      allelic_BAM_input/filtered_bam/sample4.filtered.bam.bai \
      allelic_BAM_input/filtered_bam/sample5.filtered.bam.bai \
      allelic_BAM_input/filtered_bam/sample6.filtered.bam.bai \
      allelic_BAM_input/Sambamba/sample1.markdup.txt \
      allelic_BAM_input/Sambamba/sample2.markdup.txt \
      allelic_BAM_input/Sambamba/sample3.markdup.txt \
      allelic_BAM_input/Sambamba/sample4.markdup.txt \
      allelic_BAM_input/Sambamba/sample5.markdup.txt \
      allelic_BAM_input/Sambamba/sample6.markdup.txt \
      allelic_BAM_input/bamCoverage/allele_specific/sample1.genome1.seq_depth_norm.bw \
      allelic_BAM_input/bamCoverage/allele_specific/sample2.genome1.seq_depth_norm.bw \
      allelic_BAM_input/bamCoverage/allele_specific/sample3.genome1.seq_depth_norm.bw \
      allelic_BAM_input/bamCoverage/allele_specific/sample4.genome1.seq_depth_norm.bw \
      allelic_BAM_input/bamCoverage/allele_specific/sample5.genome1.seq_depth_norm.bw \
      allelic_BAM_input/bamCoverage/allele_specific/sample6.genome1.seq_depth_norm.bw
mkdir -p output
touch /tmp/genes.gtf /tmp/genome.fa /tmp/genome.fa.fai /tmp/rmsk.txt /tmp/genes.bed /tmp/spikein_genes.gtf /tmp/genome.2bit
mkdir -p allelic_input
mkdir -p allelic_input/Ngenome
touch allelic_input/file.vcf.gz allelic_input/snpfile.txt
cp .ci_stuff/genome.fa .ci_stuff/genome.fa.fai /tmp/

# Ensure an empty snakePipes config doesn't muck anything up
snakePipes config

# createIndices
WC=`createIndices -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --genome ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtf ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz blah | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 142 ]; then exit 1 ; fi
WC=`createIndices -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --genome ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtf ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz --rmskURL http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz blah | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 148 ]; then exit 1 ; fi
WC=`createIndices -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --DAG --genome ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtf ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz --rmskURL http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz blah | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 148 ]; then exit 1 ; fi
# spikein
WC=`createIndices -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --DAG --genomeURL ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtfURL ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz --spikeinGenomeURL ftp://ftp.ensembl.org/pub/release-79/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa.gz --spikeinGtfURL  ftp://ftp.ensembl.org/pub/release-96/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.22.96.gtf.gz --rmskURL http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz blah | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 181 ]; then exit 1 ; fi


# DNA mapping
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp " | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 769 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properPairs | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 825 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properPairs --bcExtract | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 793 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --UMIDedup --properPairs --bcExtract | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 843 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --UMIDedup --properPairs | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 875 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --DAG --trim --mapq 20 --UMIDedup --properPairs | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 875 ]; then exit 1 ; fi
WC=`DNA-mapping -i SE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 693 ]; then exit 1 ; fi
WC=`DNA-mapping -i SE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properPairs | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 749 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --DAG --trim --aligner bwa | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 763 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --DAG --trim --aligner bwa-mem2 | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 763 ]; then exit 1 ; fi
#allelic
WC=`DNA-mapping -m allelic-mapping -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1,strain2 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1341 ]; then exit 1 ; fi
WC=`DNA-mapping -m allelic-mapping -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --SNPfile allelic_input/snpfile.txt --NMaskedIndex allelic_input/Ngenome .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1324 ]; then exit 1 ; fi
WC=`DNA-mapping -m allelic-mapping -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1341 ]; then exit 1 ; fi

# ChIP-seq
WC=`ChIP-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 403 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 403 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --singleEnd .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 401 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --bigWigType log2ratio .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 365 ]; then exit 1 ; fi
# fromBAM
WC=`ChIP-seq -d outdir --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 717 ]; then exit 1 ; fi
# spikein
WC=`ChIP-seq -d BAM_input --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp"  .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 667 ]; then exit 1 ; fi
# fromBAM and spikein
WC=`ChIP-seq -d outdir --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 825 ]; then exit 1 ; fi
WC=`ChIP-seq -d outdir --useSpikeInForNorm --getSizeFactorsFrom TSS --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 625 ]; then exit 1 ; fi
WC=`ChIP-seq -d outdir --useSpikeInForNorm --getSizeFactorsFrom input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 613 ]; then exit 1 ; fi
# allelic
WC=`ChIP-seq -d allelic_BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp"  .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 333 ]; then exit 1 ; fi

# mRNA-seq
WC=`mRNA-seq -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 923 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 933 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --rMats --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 949 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 643 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 989 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment-free,deepTools_qc" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1050 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --bcExtract --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 957 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --bcExtract --UMIDedup --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1007 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 848 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 567 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 904 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment-free,deepTools_qc" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 965 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --fastqc .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1016 ]; then exit 1 ; fi
WC=`mRNA-seq -i BAM_input/filtered_bam -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 657 ]; then exit 1 ; fi
#multiple comparison groups
WC=`mRNA-seq --mode alignment,alignment-free -i PE_input -o output --rMats --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 913 ]; then exit 1 ; fi
# three prime sequencing
WC=`mRNA-seq -i PE_input -o output --mode three-prime-seq --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 656 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --mode three-prime-seq,deepTools_qc --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1243 ]; then exit 1 ; fi
#allelic
WC=`mRNA-seq -m allelic-mapping,deepTools_qc -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1,strain2 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1395 ]; then exit 1 ; fi
WC=`mRNA-seq -m allelic-mapping,deepTools_qc -i allelic_BAM_input/filtered_bam --fromBAM --bamExt '.filtered.bam' -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --SNPfile allelic_input/snpfile.txt --NMaskedIndex allelic_input/Ngenome .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1118 ]; then exit 1 ; fi
WC=`mRNA-seq -m allelic-mapping,deepTools_qc -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1,strain2 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1405 ]; then exit 1 ; fi
WC=`mRNA-seq -m allelic-mapping,deepTools_qc -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --SNPfile allelic_input/snpfile.txt --NMaskedIndex allelic_input/Ngenome .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1388 ]; then exit 1 ; fi
WC=`mRNA-seq -m allelic-mapping,deepTools_qc -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1405 ]; then exit 1 ; fi
WC=`mRNA-seq -m allelic-mapping,deepTools_qc,alignment-free -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1841 ]; then exit 1 ; fi

WC=`noncoding-RNA-seq -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 755 ]; then exit 1 ; fi
WC=`noncoding-RNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 771 ]; then exit 1 ; fi
WC=`noncoding-RNA-seq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 686 ]; then exit 1 ; fi
WC=`noncoding-RNA-seq -i BAM_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 569 ]; then exit 1 ; fi
#multiple comparison groups
WC=`noncoding-RNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 793 ]; then exit 1 ; fi

# scRNA-seq
#WC=`scRNAseq -i PE_input -o output --mode Gruen --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
#if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1038 ]; then exit 1 ; fi
#WC=`scRNAseq -i PE_input -o output --mode Gruen --snakemakeOptions " --dryrun --conda-prefix /tmp" --skipRaceID --splitLib .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
#if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1015 ]; then exit 1 ; fi
WC=`scRNAseq -i PE_input -o output --mode STARsolo --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1009 ]; then exit 1 ; fi
WC=`scRNAseq -i PE_input -o output --mode STARsolo --skipVelocyto --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 901 ]; then exit 1 ; fi
WC=`scRNAseq -i PE_input -o output --mode Alevin --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 450 ]; then exit 1 ; fi
WC=`scRNAseq -i PE_input -o output --mode Alevin --skipVelocyto --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 392 ]; then exit 1 ; fi

# WGBS
WC=`WGBS -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 789 ]; then exit 1 ; fi
WC=`WGBS -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --aligner bwameth2 --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 789 ]; then exit 1 ; fi
WC=`WGBS -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --GCbias .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 798 ]; then exit 1 ; fi
WC=`WGBS -i BAM_input/filtered_bam -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --fromBAM --snakemakeOptions " --dryrun --conda-prefix /tmp" --GCbias .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 590 ]; then exit 1 ; fi
WC=`WGBS -i BAM_input/filtered_bam -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --fromBAM --fastqc --snakemakeOptions " --dryrun --conda-prefix /tmp" --GCbias .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 590 ]; then exit 1 ; fi
WC=`WGBS -i BAM_input/filtered_bam -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --fromBAM --skipBamQC --snakemakeOptions " --dryrun --conda-prefix /tmp" --GCbias .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 310 ]; then exit 1 ; fi

# ATAC-seq
WC=`ATAC-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 408 ]; then exit 1 ; fi
WC=`ATAC-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 474 ]; then exit 1 ; fi
WC=`ATAC-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller HMMRATAC .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 466 ]; then exit 1 ; fi
WC=`ATAC-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --maxFragmentSize 120 --qval 0.1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 408 ]; then exit 1 ; fi
WC=`ATAC-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 678 ]; then exit 1 ; fi

# HiC
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --correctionMethod ICE .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 533 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 501 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 557 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --enzyme DpnII .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 501 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --noTAD .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 451 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --aligner bwa-mem2 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 501 ]; then exit 1 ; fi

# preprocessing
WC=`preprocessing -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp"  --fastqc --optDedupDist 2500 | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 424 ]; then exit 1 ; fi
WC=`preprocessing -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp"  --DAG --fastqc --optDedupDist 2500 | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 424 ]; then exit 1 ; fi

rm -rf SE_input PE_input BAM_input output allelic_input allelic_BAM_input /tmp/genes.gtf /tmp/genome.fa /tmp/genome.fa.fai /tmp/rmsk.txt /tmp/genes.bed /tmp/spikein_genes.gtf
