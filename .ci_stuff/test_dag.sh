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
WC=`createIndices -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --genome ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtf ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz blah | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 229 ]; then exit 1 ; fi
WC=`createIndices -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --genome ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtf ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz --rmskURL http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz blah | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 237 ]; then exit 1 ; fi
WC=`createIndices -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --DAG --genome ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtf ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz --rmskURL http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz blah | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 237 ]; then exit 1 ; fi
# spikein
WC=`createIndices -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --DAG --genomeURL ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz --gtfURL ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz --spikeinGenomeURL ftp://ftp.ensembl.org/pub/release-79/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa.gz --spikeinGtfURL  ftp://ftp.ensembl.org/pub/release-96/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.22.96.gtf.gz --rmskURL http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz blah | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 284 ]; then exit 1 ; fi


# DNA mapping
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp " | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1420 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properPairs | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1521 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properPairs --bcExtract | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1456 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --UMIDedup --properPairs --bcExtract | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1557 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --UMIDedup --properPairs | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1622 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --DAG --trim --mapq 20 --UMIDedup --properPairs | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1622 ]; then exit 1 ; fi
WC=`DNA-mapping -i SE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1272 ]; then exit 1 ; fi
WC=`DNA-mapping -i SE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --mapq 20 --dedup --properPairs | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1373 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --DAG --trim --aligner bwa | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1402 ]; then exit 1 ; fi
WC=`DNA-mapping -i PE_input -o output .ci_stuff/organism.yaml --snakemakeOptions " --dryrun --conda-prefix /tmp" --DAG --trim --aligner bwa-mem2 | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1402 ]; then exit 1 ; fi
#allelic
WC=`DNA-mapping -m allelic-mapping -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1,strain2 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 2466 ]; then exit 1 ; fi
WC=`DNA-mapping -m allelic-mapping -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --SNPfile allelic_input/snpfile.txt --NMaskedIndex allelic_input/Ngenome .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 2445 ]; then exit 1 ; fi
WC=`DNA-mapping -m allelic-mapping -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 2466 ]; then exit 1 ; fi

# ChIP-seq
WC=`ChIP-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 630 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 609 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --singleEnd .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 628 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --bigWigType log2ratio .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 562 ]; then exit 1 ; fi
#noInput
WC=`ChIP-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 403 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/organism.yaml .ci_stuff/ChIP.sample__noControl_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
# fromBAM
WC=`ChIP-seq -d outdir --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1187 ]; then exit 1 ; fi
# fromBam and noInput
WC=`ChIP-seq -d outdir --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 801 ]; then exit 1 ; fi
# spikein
WC=`ChIP-seq -d BAM_input --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp"  .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1142 ]; then exit 1 ; fi
# spikein and noInput
WC=`ChIP-seq -d BAM_input --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp"  .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 698 ]; then exit 1 ; fi
# fromBAM and spikein
WC=`ChIP-seq -d outdir --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1447 ]; then exit 1 ; fi
WC=`ChIP-seq -d outdir --useSpikeInForNorm --getSizeFactorsFrom TSS --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1178 ]; then exit 1 ; fi
WC=`ChIP-seq -d outdir --useSpikeInForNorm --getSizeFactorsFrom input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1149 ]; then exit 1 ; fi
# allelic
WC=`ChIP-seq -d allelic_BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp"  .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_short_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 423 ]; then exit 1 ; fi
#multiComp
WC=`ChIP-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 842 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 861 ]; then exit 1 ; fi
#multiComp and fromBam
WC=`ChIP-seq -d outdir --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1399 ]; then exit 1 ; fi
WC=`ChIP-seq -d outdir --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich --fromBAM BAM_input/filtered_bam/ .ci_stuff/organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1418 ]; then exit 1 ; fi
#multiComp and spikein
WC=`ChIP-seq -d BAM_input --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp"  .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1255 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1274 ]; then exit 1 ; fi
#multiComp and spikein and noInput
WC=`ChIP-seq -d BAM_input --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp"  .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 853 ]; then exit 1 ; fi
WC=`ChIP-seq -d BAM_input --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 839 ]; then exit 1 ; fi
#multiComp and spikein and fromBam
WC=`ChIP-seq -d outdir --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1560 ]; then exit 1 ; fi
WC=`ChIP-seq -d outdir --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1045 ]; then exit 1 ; fi
#multiComp and spikein and fromBam and noInput
WC=`ChIP-seq -d outdir --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1059 ]; then exit 1 ; fi
WC=`ChIP-seq -d outdir --useSpikeInForNorm --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich --fromBAM BAM_input/filtered_bam/ .ci_stuff/spikein_organism.yaml .ci_stuff/ChIP.sample_noControl_config.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1045 ]; then exit 1 ; fi

# mRNA-seq
WC=`mRNA-seq -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1181 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1193 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --rMats --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1213 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 841 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1261 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment-free,deepTools_qc" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1330 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --bcExtract --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1217 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --bcExtract --UMIDedup --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1279 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1082 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 741 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment,deepTools_qc" --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1150 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" -m "alignment-free,deepTools_qc" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1219 ]; then exit 1 ; fi
WC=`mRNA-seq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --fastqc .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1286 ]; then exit 1 ; fi
WC=`mRNA-seq -i BAM_input/filtered_bam -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 833 ]; then exit 1 ; fi
#multiple comparison groups
WC=`mRNA-seq --mode alignment,alignment-free -i PE_input -o output --rMats --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1154 ]; then exit 1 ; fi
# three prime sequencing
WC=`mRNA-seq -i PE_input -o output --mode three-prime-seq --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1259 ]; then exit 1 ; fi
WC=`mRNA-seq -i PE_input -o output --mode three-prime-seq,deepTools_qc --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1611 ]; then exit 1 ; fi
#allelic
WC=`mRNA-seq -m allelic-mapping,deepTools_qc -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1,strain2 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1771 ]; then exit 1 ; fi
WC=`mRNA-seq -m allelic-mapping,deepTools_qc -i allelic_BAM_input/filtered_bam --fromBAM --bamExt '.filtered.bam' -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --SNPfile allelic_input/snpfile.txt --NMaskedIndex allelic_input/Ngenome .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1408 ]; then exit 1 ; fi
WC=`mRNA-seq -m allelic-mapping,deepTools_qc -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1,strain2 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1783 ]; then exit 1 ; fi
WC=`mRNA-seq -m allelic-mapping,deepTools_qc -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --SNPfile allelic_input/snpfile.txt --NMaskedIndex allelic_input/Ngenome .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1762 ]; then exit 1 ; fi
WC=`mRNA-seq -m allelic-mapping,deepTools_qc -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1783 ]; then exit 1 ; fi
WC=`mRNA-seq -m allelic-mapping,deepTools_qc,alignment-free -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --VCFfile allelic_input/file.vcf.gz --strains strain1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 2304 ]; then exit 1 ; fi

WC=`noncoding-RNA-seq -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 968 ]; then exit 1 ; fi
WC=`noncoding-RNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 988 ]; then exit 1 ; fi
WC=`noncoding-RNA-seq -i SE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 877 ]; then exit 1 ; fi
WC=`noncoding-RNA-seq -i BAM_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 714 ]; then exit 1 ; fi
#multiple comparison groups
WC=`noncoding-RNA-seq -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet_multiComp.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1016 ]; then exit 1 ; fi

# scRNA-seq
#WC=`scRNAseq -i PE_input -o output --mode Gruen --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
#if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1038 ]; then exit 1 ; fi
#WC=`scRNAseq -i PE_input -o output --mode Gruen --snakemakeOptions " --dryrun --conda-prefix /tmp" --skipRaceID --splitLib .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
#if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1015 ]; then exit 1 ; fi
WC=`scRNAseq -i PE_input -o output --mode STARsolo --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1277 ]; then exit 1 ; fi
WC=`scRNAseq -i PE_input -o output --mode STARsolo --skipVelocyto --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1143 ]; then exit 1 ; fi
WC=`scRNAseq -i PE_input -o output --mode Alevin --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 536 ]; then exit 1 ; fi
WC=`scRNAseq -i PE_input -o output --mode Alevin --skipVelocyto --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 446 ]; then exit 1 ; fi

# WGBS
WC=`WGBS -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1034 ]; then exit 1 ; fi
WC=`WGBS -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --aligner bwameth2 --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1034 ]; then exit 1 ; fi
WC=`WGBS -i PE_input -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim --GCbias .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 1045 ]; then exit 1 ; fi
WC=`WGBS -i BAM_input/filtered_bam -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --fromBAM --snakemakeOptions " --dryrun --conda-prefix /tmp" --GCbias .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 759 ]; then exit 1 ; fi
WC=`WGBS -i BAM_input/filtered_bam -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --fromBAM --fastqc --snakemakeOptions " --dryrun --conda-prefix /tmp" --GCbias .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 759 ]; then exit 1 ; fi
WC=`WGBS -i BAM_input/filtered_bam -o output --sampleSheet .ci_stuff/test_sampleSheet.tsv --fromBAM --skipBamQC --snakemakeOptions " --dryrun --conda-prefix /tmp" --GCbias .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 400 ]; then exit 1 ; fi

# ATAC-seq
WC=`ATAC-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 527 ]; then exit 1 ; fi
WC=`ATAC-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller Genrich .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 615 ]; then exit 1 ; fi
WC=`ATAC-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --peakCaller HMMRATAC .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 600 ]; then exit 1 ; fi
WC=`ATAC-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --maxFragmentSize 120 --qval 0.1 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 527 ]; then exit 1 ; fi
WC=`ATAC-seq -d BAM_input --sampleSheet .ci_stuff/test_sampleSheet.tsv --snakemakeOptions " --dryrun --conda-prefix /tmp" --fromBAM BAM_input/filtered_bam/ .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 869 ]; then exit 1 ; fi

# HiC
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --correctionMethod ICE .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 697 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 653 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --trim .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 721 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --enzyme DpnII .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 653 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --noTAD .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 591 ]; then exit 1 ; fi
WC=`HiC -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp" --aligner bwa-mem2 .ci_stuff/organism.yaml | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 653 ]; then exit 1 ; fi

# preprocessing
WC=`preprocessing -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp"  --fastqc --optDedupDist 2500 | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 558 ]; then exit 1 ; fi
WC=`preprocessing -i PE_input -o output --snakemakeOptions " --dryrun --conda-prefix /tmp"  --DAG --fastqc --optDedupDist 2500 | tee >(cat 1>&2) | grep -v "conda installation\|Conda environment" | sed '/^\s*$/d' | wc -l`
if [ ${PIPESTATUS[0]} -ne 0 ] || [ $WC -ne 558 ]; then exit 1 ; fi

rm -rf SE_input PE_input BAM_input output allelic_input allelic_BAM_input /tmp/genes.gtf /tmp/genome.fa /tmp/genome.fa.fai /tmp/rmsk.txt /tmp/genes.bed /tmp/spikein_genes.gtf
