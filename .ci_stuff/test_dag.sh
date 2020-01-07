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
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 886 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properPairs | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 954 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properPairs --bcExtract | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 914 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --UMIDedup --properPairs --bcExtract | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 970 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --UMIDedup --properPairs | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1010 ]; then exit 1 ; fi
WC=`DNA-mapping -i SE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 786 ]; then exit 1 ; fi
WC=`DNA-mapping -i SE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properPairs | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 854 ]; then exit 1 ; fi

# ChIP-seq
WC=`ChIP-seq -d BAM_input --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 270 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --snakemakeOptions " --dryrun --conda-prefix /tmp" --singleEnd .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 270 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --snakemakeOptions " --dryrun --conda-prefix /tmp" --bigWigType log2ratio .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 224 ]; then exit 1 ; fi

# ATAC-seq
WC=`ATAC-seq -d BAM_input --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 337 ]; then exit 1 ; fi
WC=`ATAC-seq -d BAM_input --snakemakeOptions " --dryrun --conda-prefix /tmp" --maxFragmentSize 120 --qval 0.1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 337 ]; then exit 1 ; fi

# mRNA-seq
WC=`mRNA-seq -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1006 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 734 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1074 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment-free,deepTools_qc" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1156 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --bcExtract --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1032 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --bcExtract --UMIDedup --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1088 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 895 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 634 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 963 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment-free,deepTools_qc" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1045 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --fastqc .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1099 ]; then exit 1 ; fi

# HiC
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --correctionMethod ICE .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 624 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 580 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 648 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --enzyme DpnII .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 580 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --noTAD .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 518 ]; then exit 1 ; fi

# scRNA-seq
WC=`scRNAseq -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1282 ]; then exit 1 ; fi
WC=`scRNAseq -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --skipRaceID --splitLib .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1253 ]; then exit 1 ; fi
WC=`scRNAseq -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --mode STARsolo .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1014 ]; then exit 1 ; fi

# WGBS
WC=`WGBS -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 849 ]; then exit 1 ; fi
WC=`WGBS -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --GCbias .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 928 ]; then exit 1 ; fi

# preprocessing
WC=`preprocessing -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp"  --fastqc --optDedupDist 2500 | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 488 ]; then exit 1 ; fi

# createIndices
WC=`preprocessing -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --genome ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtf ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz blah | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 488 ]; then exit 1 ; fi
WC=`preprocessing -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --genome ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtf ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz --rmskURL http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz blah | tee >(cat 1>&2) | grep -v "Conda environment" | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 488 ]; then exit 1 ; fi

rm -rf SE_input PE_input BAM_input output /tmp/genes.gtf /tmp/genome.fa /tmp/genome.fa.fai
