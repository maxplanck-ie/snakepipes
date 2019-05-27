#!/bin/bash
set -ex
if [[ ${CI:-"false"} == "true" ]]; then
    export PATH="$HOME/miniconda/bin:$PATH"
    hash -r
    python -m pip install --no-deps --ignore-installed .
fi

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
touch /tmp/genes.gtf /tmp/genome.fa /tmp/genome.fa.fai

# DNA mapping
WC=`DNA-mapping -i PE_input -o output --tempdir /tmp .ci_stuff/organism.yaml --snakemake_options " --dryrun --conda-prefix /tmp" | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 884 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output --tempdir /tmp .ci_stuff/organism.yaml --snakemake_options " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properpairs | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 952 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output --tempdir /tmp .ci_stuff/organism.yaml --snakemake_options " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properpairs --umibc | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 912 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output --tempdir /tmp .ci_stuff/organism.yaml --snakemake_options " --dryrun --conda-prefix /tmp" --trim --mapq 20 --umidedup --properpairs --umibc | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 968 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output --tempdir /tmp .ci_stuff/organism.yaml --snakemake_options " --dryrun --conda-prefix /tmp" --trim --mapq 20 --umidedup --properpairs | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1008 ]; then exit 1 ; fi
WC=`DNA-mapping -i SE_input -o output --tempdir /tmp .ci_stuff/organism.yaml --snakemake_options " --dryrun --conda-prefix /tmp" | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 784 ]; then exit 1 ; fi
WC=`DNA-mapping -i SE_input -o output --tempdir /tmp .ci_stuff/organism.yaml --snakemake_options " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properpairs | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 852 ]; then exit 1 ; fi

# ChIP-seq
WC=`ChIP-seq -d BAM_input --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 268 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" --single-end .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 268 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" --bigWigType log2ratio .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 222 ]; then exit 1 ; fi

# ATAC-seq
WC=`ATAC-seq -d BAM_input --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 336 ]; then exit 1 ; fi
WC=`ATAC-seq -d BAM_input --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" --fragmentSize-cutoff 120 --qval 0.1 .ci_stuff/organism.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 336 ]; then exit 1 ; fi

# RNA-seq
WC=`RNA-seq -i PE_input -o output --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 989 ]; then exit 1 ; fi
WC=`RNA-seq -i PE_input -o output --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" -m "alignment" .ci_stuff/organism.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 678 ]; then exit 1 ; fi
WC=`RNA-seq -i PE_input -o output --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --trim .ci_stuff/organism.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1018 ]; then exit 1 ; fi
WC=`RNA-seq -i PE_input -o output --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --umibc --trim .ci_stuff/organism.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 978 ]; then exit 1 ; fi #TEST?!!
WC=`RNA-seq -i SE_input -o output --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 878 ]; then exit 1 ; fi
WC=`RNA-seq -i SE_input -o output --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" -m "alignment" .ci_stuff/organism.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 578 ]; then exit 1 ; fi
WC=`RNA-seq -i SE_input -o output --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --trim .ci_stuff/organism.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 907 ]; then exit 1 ; fi

# HiC - fails
WC=`HiC -i PE_input -o output --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 622 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" --trim .ci_stuff/organism.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 690 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" --enzyme DpnII .ci_stuff/organism.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 622 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" --noTAD .ci_stuff/organism.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 560 ]; then exit 1 ; fi

# scRNA-seq
WC=`scRNAseq -i PE_input -o output --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1180 ]; then exit 1 ; fi
WC=`scRNAseq -i PE_input -o output --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" --skipRaceID --split_lib .ci_stuff/organism.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1151 ]; then exit 1 ; fi

# WGBS
WC=`WGBS -i PE_input -o output --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1061 ]; then exit 1 ; fi
WC=`WGBS -i PE_input -o output --tempdir /tmp --snakemake_options " --dryrun --conda-prefix /tmp" --skipDOC --skipGCbias --nextera .ci_stuff/organism.yaml | tee /dev/stderr | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 943 ]; then exit 1 ; fi

rm -rf SE_input PE_input BAM_input output /tmp/genes.gtf /tmp/genome.fa /tmp/genome.fa.fai
