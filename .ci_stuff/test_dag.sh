#!/bin/bash
set -ex
export PATH="$HOME/miniconda/bin:$PATH"
hash -r
python -m pip install --no-deps --ignore-installed .

# Needed by DNA, HiC, RNA-seq, WGBS and scRNA-seq workflows
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
mkdir -p BAM_input/deepTools_qc/bamPEFragmentSize BAM_input/filtered_bam BAM_input/Sambamba
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
      BAM_input/deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv
mkdir -p output

# DNA mapping
WC=`DNA-mapping -i PE_input -o output --tempdir /tmp mm10 --snakemake_options " --dryrun" | tee /dev/stderr | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 734 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output --tempdir /tmp mm10 --snakemake_options " --dryrun" --trim --mapq 20 --dedup --properpairs | tee /dev/stderr | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 936 ]; then exit 1 ; fi
WC=`DNA-mapping -i SE_input -o output --tempdir /tmp mm10 --snakemake_options " --dryrun" | tee /dev/stderr | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 684 ]; then exit 1 ; fi
WC=`DNA-mapping -i SE_input -o output --tempdir /tmp mm10 --snakemake_options " --dryrun" --trim --mapq 20 --dedup --properpairs | tee /dev/stderr | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 820 ]; then exit 1 ; fi

# ChIP-seq
WC=`ChIP-seq -d BAM_input --tempdir /tmp --snakemake_options " --dryrun" mm10 .ci_stuff/ChIP.sample_config.yaml | tee /dev/stderr | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 268 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --tempdir /tmp --snakemake_options " --dryrun" --single-end mm10 .ci_stuff/ChIP.sample_config.yaml | tee /dev/stderr | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 268 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --tempdir /tmp --snakemake_options " --dryrun" --bigWigType log2ratio mm10 .ci_stuff/ChIP.sample_config.yaml | tee /dev/stderr | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 222 ]; then exit 1 ; fi

rm -rf SE_input PE_input BAM_input output
