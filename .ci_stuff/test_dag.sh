#!/bin/bash
set -ex

# Needed by DNA, HiC, mRNAseq, WGBS and scRNAseq workflows
mkdir -p PE_input
touch PE_input/sample1_R1.fastq.gz PE_input/sample1_R2.fastq.gz \
      PE_input/sample2_R1.fastq.gz PE_input/sample2_R2.fastq.gz \
      PE_input/sample3_R1.fastq.gz PE_input/sample3_R2.fastq.gz \
      PE_input/sample4_R1.fastq.gz PE_input/sample4_R2.fastq.gz \
      PE_input/sample5_R1.fastq.gz PE_input/sample5_R2.fastq.gz \
      PE_input/sample6_R1.fastq.gz PE_input/sample6_R2.fastq.gz \
      PE_input/sample7_R1.fastq.gz PE_input/sample7_R2.fastq.gz \
      PE_input/sample8_R1.fastq.gz PE_input/sample8_R2.fastq.gz \
      PE_input/sample9_R1.fastq.gz PE_input/sample9_R2.fastq.gz
mkdir -p SE_input
touch SE_input/sample1_R1.fastq.gz \
      SE_input/sample2_R1.fastq.gz \
      SE_input/sample3_R1.fastq.gz \
      SE_input/sample4_R1.fastq.gz \
      SE_input/sample5_R1.fastq.gz \
      SE_input/sample6_R1.fastq.gz \
      SE_input/sample7_R1.fastq.gz \
      SE_input/sample8_R1.fastq.gz \
      SE_input/sample9_R1.fastq.gz
# Needed by ChIP and ATAC workflows
mkdir -p BAM_input/deepTools_qc/bamPEFragmentSize BAM_input/filtered_bam BAM_input/Sambamba BAM_input/bamCoverage
touch BAM_input/sample1.bam \
      BAM_input/sample2.bam \
      BAM_input/sample3.bam \
      BAM_input/sample4.bam \
      BAM_input/sample5.bam \
      BAM_input/sample6.bam \
      BAM_input/sample7.bam \
      BAM_input/sample8.bam \
      BAM_input/sample9.bam \
      BAM_input/filtered_bam/sample1.filtered.bam \
      BAM_input/filtered_bam/sample2.filtered.bam \
      BAM_input/filtered_bam/sample3.filtered.bam \
      BAM_input/filtered_bam/sample4.filtered.bam \
      BAM_input/filtered_bam/sample5.filtered.bam \
      BAM_input/filtered_bam/sample6.filtered.bam \
      BAM_input/filtered_bam/sample7.filtered.bam \
      BAM_input/filtered_bam/sample8.filtered.bam \
      BAM_input/filtered_bam/sample9.filtered.bam \
      BAM_input/filtered_bam/sample1.filtered.bam.bai \
      BAM_input/filtered_bam/sample2.filtered.bam.bai \
      BAM_input/filtered_bam/sample3.filtered.bam.bai \
      BAM_input/filtered_bam/sample4.filtered.bam.bai \
      BAM_input/filtered_bam/sample5.filtered.bam.bai \
      BAM_input/filtered_bam/sample6.filtered.bam.bai \
      BAM_input/filtered_bam/sample7.filtered.bam.bai \
      BAM_input/filtered_bam/sample8.filtered.bam.bai \
      BAM_input/filtered_bam/sample9.filtered.bam.bai \
      BAM_input/Sambamba/sample1.markdup.txt \
      BAM_input/Sambamba/sample2.markdup.txt \
      BAM_input/Sambamba/sample3.markdup.txt \
      BAM_input/Sambamba/sample4.markdup.txt \
      BAM_input/Sambamba/sample5.markdup.txt \
      BAM_input/Sambamba/sample6.markdup.txt \
      BAM_input/Sambamba/sample7.markdup.txt \
      BAM_input/Sambamba/sample8.markdup.txt \
      BAM_input/Sambamba/sample9.markdup.txt \
      BAM_input/deepTools_qc/bamPEFragmentSize/fragmentSize.metric.tsv \
      BAM_input/bamCoverage/sample1.filtered.seq_depth_norm.bw \
      BAM_input/bamCoverage/sample2.filtered.seq_depth_norm.bw \
      BAM_input/bamCoverage/sample3.filtered.seq_depth_norm.bw \
      BAM_input/bamCoverage/sample4.filtered.seq_depth_norm.bw \
      BAM_input/bamCoverage/sample5.filtered.seq_depth_norm.bw \
      BAM_input/bamCoverage/sample6.filtered.seq_depth_norm.bw \
      BAM_input/bamCoverage/sample7.filtered.seq_depth_norm.bw \
      BAM_input/bamCoverage/sample8.filtered.seq_depth_norm.bw \
      BAM_input/bamCoverage/sample9.filtered.seq_depth_norm.bw 

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
      allelic_BAM_input/allelic_bams/sample1.allele_flagged.sorted.bam \
      allelic_BAM_input/allelic_bams/sample1.unassigned.sorted.bam \
      allelic_BAM_input/allelic_bams/sample2.allele_flagged.sorted.bam \
      allelic_BAM_input/allelic_bams/sample2.unassigned.sorted.bam \
      allelic_BAM_input/allelic_bams/sample3.allele_flagged.sorted.bam \
      allelic_BAM_input/allelic_bams/sample3.unassigned.sorted.bam \
      allelic_BAM_input/allelic_bams/sample4.allele_flagged.sorted.bam \
      allelic_BAM_input/allelic_bams/sample4.unassigned.sorted.bam \
      allelic_BAM_input/allelic_bams/sample5.allele_flagged.sorted.bam \
      allelic_BAM_input/allelic_bams/sample5.unassigned.sorted.bam \
      allelic_BAM_input/allelic_bams/sample6.allele_flagged.sorted.bam \
      allelic_BAM_input/allelic_bams/sample6.unassigned.sorted.bam \
      allelic_BAM_input/allelic_bams/sample1.allele_flagged.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample1.unassigned.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample2.allele_flagged.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample2.unassigned.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample3.allele_flagged.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample3.unassigned.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample4.allele_flagged.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample4.unassigned.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample5.allele_flagged.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample5.unassigned.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample6.allele_flagged.sorted.bam.bai \
      allelic_BAM_input/allelic_bams/sample6.unassigned.sorted.bam.bai \
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
mkdir -p /tmp/SalmonIndex /tmp/annotation
touch /tmp/SalmonIndex/decoys.txt
touch /tmp/annotation/cDNA_introns.joint.t2g
touch /tmp/genes.t2g

# Ensure an empty snakePipes config doesn't muck anything up
snakePipes config --tempDir /tmp

# createIndices
WC=`createIndices -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --genome ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtf ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz blah | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 213 ]; then exit 1 ; fi
WC=`createIndices -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --genome ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtf ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz --rmskURL http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz blah | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 221 ]; then exit 1 ; fi
WC=`createIndices -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --DAG --genome ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtf ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz --rmskURL http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz blah | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 221 ]; then exit 1 ; fi
# spikein
WC=`createIndices -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --DAG --genomeURL ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtfURL ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz --spikeinGenomeURL ftp://ftp.ensembl.org/pub/release-79/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa.gz --spikeinGtfURL  ftp://ftp.ensembl.org/pub/release-96/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.22.96.gtf.gz --rmskURL http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz blah | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 268 ]; then exit 1 ; fi


# DNA mapping
WC=`DNAmapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp " | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1332 ]; then exit 1 ; fi
WC=`DNAmapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properPairs | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1424 ]; then exit 1 ; fi
WC=`DNAmapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properPairs --bcExtract | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1350 ]; then exit 1 ; fi
WC=`DNAmapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --UMIDedup --properPairs --bcExtract | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1433 ]; then exit 1 ; fi
WC=`DNAmapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --UMIDedup --properPairs | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1507 ]; then exit 1 ; fi
WC=`DNAmapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --DAG --trim --mapq 20 --UMIDedup --properPairs | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1507 ]; then exit 1 ; fi
WC=`DNAmapping -i SE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1184 ]; then exit 1 ; fi
WC=`DNAmapping -i SE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properPairs | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1276 ]; then exit 1 ; fi
WC=`DNAmapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --DAG --trim --aligner bwa | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1314 ]; then exit 1 ; fi
WC=`DNAmapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --DAG --trim --aligner bwa-mem2 | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1314 ]; then exit 1 ; fi
#allelic
WC=`DNAmapping -m allelic-mapping -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1,strain2 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 2273 ]; then exit 1 ; fi
WC=`DNAmapping -m allelic-mapping -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --SNPfile allelic_input/snpfile.txt --NMaskedIndex allelic_input/Ngenome .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 2254 ]; then exit 1 ; fi
WC=`DNAmapping -m allelic-mapping -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 2273 ]; then exit 1 ; fi

# ChIPseq
WC=`ChIPseq -d BAM_input --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 377 ]; then exit 1 ; fi
WC=`ChIPseq -d BAM_input --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 365 ]; then exit 1 ; fi
WC=`ChIPseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 580 ]; then exit 1 ; fi
WC=`ChIPseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_broad_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 778 ]; then exit 1 ; fi
WC=`ChIPseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 556 ]; then exit 1 ; fi
WC=`ChIPseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --singleEnd .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 578 ]; then exit 1 ; fi
WC=`ChIPseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --bigWigType log2ratio .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 518 ]; then exit 1 ; fi
WC=`ChIPseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller SEACR .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 853 ]; then exit 1 ; fi
WC=`ChIPseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller SEACR --useSpikeInForNorm .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1256 ]; then exit 1 ; fi
#noInput
WC=`ChIPseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 374 ]; then exit 1 ; fi
WC=`ChIPseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 320 ]; then exit 1 ; fi
WC=`ChIPseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller SEACR .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 590 ]; then exit 1 ; fi
WC=`ChIPseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller SEACR --useSpikeInForNorm .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 796 ]; then exit 1 ; fi
# fromBAM
WC=`ChIPseq -d outdir --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1104 ]; then exit 1 ; fi
# fromBam and noInput
WC=`ChIPseq -d outdir --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 748 ]; then exit 1 ; fi
# spikein
WC=`ChIPseq -d BAM_input --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp"  .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1045 ]; then exit 1 ; fi
# spikein and noInput
WC=`ChIPseq -d BAM_input --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp"  .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 642 ]; then exit 1 ; fi
# fromBAM and spikein
WC=`ChIPseq -d outdir --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1341 ]; then exit 1 ; fi
WC=`ChIPseq -d outdir --useSpikeInForNorm --getSizeFactorsFrom TSS --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1098 ]; then exit 1 ; fi
WC=`ChIPseq -d outdir --useSpikeInForNorm --getSizeFactorsFrom input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1070 ]; then exit 1 ; fi
# allelic
WC=`ChIPseq -d allelic_BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp"  .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_short_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 391 ]; then exit 1 ; fi
#multiComp
WC=`ChIPseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 772 ]; then exit 1 ; fi
WC=`ChIPseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 784 ]; then exit 1 ; fi
WC=`ChIPseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_broad_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 970 ]; then exit 1 ; fi
#multiComp and fromBam
WC=`ChIPseq -d outdir --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1296 ]; then exit 1 ; fi
WC=`ChIPseq -d outdir --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich --fromBAM BAM_input/filtered_bam/ .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1308 ]; then exit 1 ; fi
#multiComp and spikein
WC=`ChIPseq -d BAM_input --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp"  .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1220 ]; then exit 1 ; fi
WC=`ChIPseq -d BAM_input --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1232 ]; then exit 1 ; fi
#multiComp and spikein and noInput
WC=`ChIPseq -d BAM_input --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp"  .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 783 ]; then exit 1 ; fi
WC=`ChIPseq -d BAM_input --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 765 ]; then exit 1 ; fi
#multiComp and spikein and fromBam
WC=`ChIPseq -d outdir --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1516 ]; then exit 1 ; fi
WC=`ChIPseq -d outdir --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 965 ]; then exit 1 ; fi
#multiComp and spikein and fromBam and noInput
WC=`ChIPseq -d outdir --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 983 ]; then exit 1 ; fi
WC=`ChIPseq -d outdir --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 965 ]; then exit 1 ; fi
#multiComp and noInput and fromBam
WC=`ChIPseq -d outdir --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller SEACR --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1105 ]; then exit 1 ; fi


# mRNAseq
WC=`mRNAseq -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1563 ]; then exit 1 ; fi
WC=`mRNAseq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1574 ]; then exit 1 ; fi
WC=`mRNAseq -i PE_input -o output --rMats --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1593 ]; then exit 1 ; fi
WC=`mRNAseq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1133 ]; then exit 1 ; fi
WC=`mRNAseq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1666 ]; then exit 1 ; fi
WC=`mRNAseq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment-free,deepTools_qc" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1754 ]; then exit 1 ; fi
WC=`mRNAseq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --bcExtract --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1592 ]; then exit 1 ; fi
WC=`mRNAseq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --bcExtract --UMIDedup --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1666 ]; then exit 1 ; fi
WC=`mRNAseq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1416 ]; then exit 1 ; fi
WC=`mRNAseq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 985 ]; then exit 1 ; fi
WC=`mRNAseq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1508 ]; then exit 1 ; fi
WC=`mRNAseq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment-free,deepTools_qc" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1596 ]; then exit 1 ; fi
WC=`mRNAseq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --fastqc .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1692 ]; then exit 1 ; fi
WC=`mRNAseq -i BAM_input/filtered_bam -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1067 ]; then exit 1 ; fi
#multiple comparison groups
WC=`mRNAseq --mode alignment,alignment-free -i PE_input -o output --rMats --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1503 ]; then exit 1 ; fi
# three prime sequencing
WC=`mRNAseq -i PE_input -o output --mode three-prime-seq --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1692 ]; then exit 1 ; fi
WC=`mRNAseq -i PE_input -o output --mode three-prime-seq,deepTools_qc --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 2133 ]; then exit 1 ; fi
#allelic
WC=`mRNAseq -m allelic-mapping,deepTools_qc -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1,strain2 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 2331 ]; then exit 1 ; fi
WC=`mRNAseq -m allelic-mapping,deepTools_qc -i allelic_BAM_input/filtered_bam --fromBAM --bamExt '.filtered.bam' -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --SNPfile allelic_input/snpfile.txt --NMaskedIndex allelic_input/Ngenome .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1291 ]; then exit 1 ; fi
WC=`mRNAseq -m allelic-mapping,deepTools_qc -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1,strain2 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 2342 ]; then exit 1 ; fi
WC=`mRNAseq -m allelic-mapping,deepTools_qc -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --SNPfile allelic_input/snpfile.txt --NMaskedIndex allelic_input/Ngenome .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 2323 ]; then exit 1 ; fi
WC=`mRNAseq -m allelic-mapping,deepTools_qc -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 2342 ]; then exit 1 ; fi
WC=`mRNAseq -m allelic-mapping,deepTools_qc,alignment-free -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 3039 ]; then exit 1 ; fi
WC=`mRNAseq -m allelic-counting -i allelic_BAM_input/allelic_bams --fromBAM --bamExt '.sorted.bam' -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp"  .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 638 ]; then exit 1 ; fi
WC=`mRNAseq -m allelic-counting -i allelic_BAM_input/allelic_bams --fromBAM --bamExt '.sorted.bam' -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp"  .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 638 ]; then exit 1 ; fi
#allelic+multicomp
WC=`mRNAseq -m allelic-counting -i allelic_BAM_input/allelic_bams --fromBAM --bamExt '.sorted.bam' -o output --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp"  .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 668 ]; then exit 1 ; fi
WC=`mRNAseq -m allelic-mapping,deepTools_qc -i allelic_BAM_input/filtered_bam --fromBAM --bamExt '.filtered.bam' -o output --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --SNPfile allelic_input/snpfile.txt --NMaskedIndex allelic_input/Ngenome .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1321 ]; then exit 1 ; fi
WC=`mRNAseq -m allelic-mapping,deepTools_qc,alignment-free -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 3048 ]; then exit 1 ; fi

#ncRNAseq
WC=`ncRNAseq -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1290 ]; then exit 1 ; fi
WC=`ncRNAseq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1310 ]; then exit 1 ; fi
WC=`ncRNAseq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1152 ]; then exit 1 ; fi
WC=`ncRNAseq -i BAM_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 904 ]; then exit 1 ; fi
#multiple comparison groups
WC=`ncRNAseq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1338 ]; then exit 1 ; fi

# scRNAseq
#WC=`scRNAseq -i PE_input -o output --mode Gruen --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
#if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1038 ]; then exit 1 ; fi
#WC=`scRNAseq -i PE_input -o output --mode Gruen --snakemakeOptions " --dryrun --conda-prefix /tmp" --skipRaceID --splitLib .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
#if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1015 ]; then exit 1 ; fi
WC=`scRNAseq -i PE_input -o output --mode STARsolo --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1642 ]; then exit 1 ; fi
WC=`scRNAseq -i PE_input -o output --mode STARsolo --skipVelocyto --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1485 ]; then exit 1 ; fi
WC=`scRNAseq -i PE_input -o output --mode Alevin --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 723 ]; then exit 1 ; fi
WC=`scRNAseq -i PE_input -o output --mode Alevin --skipVelocyto --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 613 ]; then exit 1 ; fi

# WGBS
WC=`WGBS -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1337 ]; then exit 1 ; fi
WC=`WGBS -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1380 ]; then exit 1 ; fi
WC=`WGBS -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --aligner bwameth2 --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1380 ]; then exit 1 ; fi
WC=`WGBS -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --GCbias .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1390 ]; then exit 1 ; fi
WC=`WGBS -i BAM_input/filtered_bam -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --fromBAM --snakemakeOptions " --dryrun --conda-prefix /tmp" --GCbias .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 984 ]; then exit 1 ; fi
WC=`WGBS -i BAM_input/filtered_bam -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --fromBAM --fastqc --snakemakeOptions " --dryrun --conda-prefix /tmp" --GCbias .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 984 ]; then exit 1 ; fi
WC=`WGBS -i BAM_input/filtered_bam -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --fromBAM --skipBamQC --snakemakeOptions " --dryrun --conda-prefix /tmp" --GCbias .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 530 ]; then exit 1 ; fi

# ATACseq
WC=`ATACseq -d BAM_input --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 478 ]; then exit 1 ; fi
WC=`ATACseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 626 ]; then exit 1 ; fi
WC=`ATACseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 736 ]; then exit 1 ; fi
WC=`ATACseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller HMMRATAC .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 719 ]; then exit 1 ; fi
WC=`ATACseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --maxFragmentSize 120 --qval 0.1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 626 ]; then exit 1 ; fi
WC=`ATACseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1076 ]; then exit 1 ; fi
#multicomp
WC=`ATACseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 767 ]; then exit 1 ; fi
WC=`ATACseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 893 ]; then exit 1 ; fi
WC=`ATACseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller HMMRATAC .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 860 ]; then exit 1 ; fi
WC=`ATACseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --maxFragmentSize 120 --qval 0.1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 767 ]; then exit 1 ; fi
WC=`ATACseq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1217 ]; then exit 1 ; fi


# HiC
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --correctionMethod ICE .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 974 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 900 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 992 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --enzyme DpnII .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 900 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --noTAD .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 817 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --aligner bwa-mem2 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 900 ]; then exit 1 ; fi

# preprocessing
WC=`preprocessing -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp"  --fastqc --optDedupDist 2500 | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 779 ]; then exit 1 ; fi
WC=`preprocessing -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp"  --DAG --fastqc --optDedupDist 2500 | tee >(cat 1>&2) | grep -v -f .ci_stuff/test_ignore_patterns.txt | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 779 ]; then exit 1 ; fi

rm -rf SE_input PE_input BAM_input output allelic_input allelic_BAM_input /tmp/genes.gtf /tmp/genome.fa /tmp/genome.fa.fai /tmp/rmsk.txt /tmp/genes.bed /tmp/spikein_genes.gtf
