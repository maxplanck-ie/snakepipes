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
touch /tmp/genes.gtf /tmp/genome.fa /tmp/genome.fa.fai /tmp/rmsk.txt /tmp/genes.bed
mkdir -p allelic_input
mkdir -p allelic_input/Ngenome
touch allelic_input/file.vcf.gz allelic_input/snpfile.txt

# Ensure an empty snakePipes config doesn't muck anything up
snakePipes config

# DNA mapping
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp " | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 701 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properPairs | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 757 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properPairs --bcExtract | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 728 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --UMIDedup --properPairs --bcExtract | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 772 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --UMIDedup --properPairs | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 801 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --DAG --trim --mapq 20 --UMIDedup --properPairs | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 801 ]; then exit 1 ; fi
WC=`DNA-mapping -i SE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 625 ]; then exit 1 ; fi
WC=`DNA-mapping -i SE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properPairs | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 681 ]; then exit 1 ; fi
#allelic
WC=`DNA-mapping -m allelic-mapping -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1,strain2 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1365 ]; then exit 1 ; fi
WC=`DNA-mapping -m allelic-mapping -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --SNPfile allelic_input/snpfile.txt --NMaskedIndex allelic_input/Ngenome .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1348 ]; then exit 1 ; fi
WC=`DNA-mapping -m allelic-mapping -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1365 ]; then exit 1 ; fi

# ChIP-seq
WC=`ChIP-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 431 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 381 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --singleEnd .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 431 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --bigWigType log2ratio .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 393 ]; then exit 1 ; fi
# fromBAM
WC=`ChIP-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 701 ]; then exit 1 ; fi



# mRNA-seq
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 799 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 575 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 855 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment-free,deepTools_qc" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 907 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --bcExtract --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 825 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --bcExtract --UMIDedup --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 869 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 714 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 499 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 770 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment-free,deepTools_qc" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 822 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --fastqc .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 882 ]; then exit 1 ; fi
WC=`mRNA-seq -i BAM_input/filtered_bam -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 578 ]; then exit 1 ; fi
#allelic
WC=`mRNA-seq -m allelic-mapping -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1,strain2 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 992 ]; then exit 1 ; fi
WC=`mRNA-seq -m allelic-mapping -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --SNPfile allelic_input/snpfile.txt --NMaskedIndex allelic_input/Ngenome .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 975 ]; then exit 1 ; fi
WC=`mRNA-seq -m allelic-mapping -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 992 ]; then exit 1 ; fi

# noncoding-RNA-seq
WC=`noncoding-RNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 655 ]; then exit 1 ; fi
WC=`noncoding-RNA-seq -i BAM_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 491 ]; then exit 1 ; fi


# scRNA-seq
WC=`scRNAseq -i PE_input -o output --mode Gruen --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1006 ]; then exit 1 ; fi
WC=`scRNAseq -i PE_input -o output --mode Gruen --snakemakeOptions " --dryrun --conda-prefix /tmp" --skipRaceID --splitLib .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1243 ]; then exit 1 ; fi
WC=`scRNAseq -i PE_input -o output --mode STARsolo --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1073 ]; then exit 1 ; fi

# WGBS
WC=`WGBS -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 891 ]; then exit 1 ; fi
WC=`WGBS -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --GCbias .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 970 ]; then exit 1 ; fi
WC=`WGBS -i BAM_input/filtered_bam -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --fromBAM --snakemakeOptions " --dryrun --conda-prefix /tmp" --GCbias .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 690 ]; then exit 1 ; fi
WC=`WGBS -i BAM_input/filtered_bam -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --fromBAM --fastqc --snakemakeOptions " --dryrun --conda-prefix /tmp" --GCbias .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 690 ]; then exit 1 ; fi

# ATAC-seq
WC=`ATAC-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 536 ]; then exit 1 ; fi
WC=`ATAC-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 554 ]; then exit 1 ; fi
WC=`ATAC-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller HMMRATAC .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 601 ]; then exit 1 ; fi
WC=`ATAC-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --maxFragmentSize 120 --qval 0.1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 536 ]; then exit 1 ; fi
WC=`ATAC-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 832 ]; then exit 1 ; fi

# HiC
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --correctionMethod ICE .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 614 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 570 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 638 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --enzyme DpnII .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 570 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --noTAD .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 508 ]; then exit 1 ; fi

# preprocessing
WC=`preprocessing -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp"  --fastqc --optDedupDist 2500 | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 488 ]; then exit 1 ; fi
WC=`preprocessing -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp"  --DAG --fastqc --optDedupDist 2500 | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 488 ]; then exit 1 ; fi

# createIndices
WC=`createIndices -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --genome ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtf ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz blah | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 143 ]; then exit 1 ; fi
WC=`createIndices -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --genome ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtf ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz --rmskURL http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz blah | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 150 ]; then exit 1 ; fi
WC=`createIndices -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --DAG --genome ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtf ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz --rmskURL http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz blah | tee >(cat 1>&2) | grep -v "Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 150 ]; then exit 1 ; fi

rm -rf SE_input PE_input BAM_input output allelic_input /tmp/genes.gtf /tmp/genome.fa /tmp/genome.fa.fai /tmp/rmsk.txt /tmp/genes.bed
