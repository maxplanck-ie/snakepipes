#!/bin/bash
export PATH="$HOME/miniconda/bin:$PATH"
hash -r
python -m pip install --no-deps --ignore-installed .


# Make some fake data
mkdir PE_input
touch PE_input/sample1_R1.fastq.gz PE_input/sample1_R2.fastq.gz \
      PE_input/sample2_R1.fastq.gz PE_input/sample2_R2.fastq.gz \
      PE_input/sample3_R1.fastq.gz PE_input/sample3_R2.fastq.gz \
      PE_input/sample4_R1.fastq.gz PE_input/sample4_R2.fastq.gz \
      PE_input/sample5_R1.fastq.gz PE_input/sample5_R2.fastq.gz \
      PE_input/sample6_R1.fastq.gz PE_input/sample6_R2.fastq.gz
mkdir SE_input
touch SE_input/sample1_R1.fastq.gz \
      SE_input/sample2_R1.fastq.gz \
      SE_input/sample3_R1.fastq.gz \
      SE_input/sample4_R1.fastq.gz \
      SE_input/sample5_R1.fastq.gz \
      SE_input/sample6_R1.fastq.gz

# DNA mapping
rm -rf output && mkdir output
DNA-mapping -i PE_input -o output mm10 --snakemake_options " --dryrun"
rm -rf output && mkdir output
DNA-mapping -i PE_input -o output mm10 --snakemake_options " --dryrun" --trim --mapq 20 --dedup --properpairs
rm -rf output && mkdir output
DNA-mapping -i SE_input -o output mm10 --snakemake_options " --dryrun"
rm -rf output && mkdir output
DNA-mapping -i SE_input -o output mm10 --snakemake_options " --dryrun" --trim --mapq 20 --dedup --properpairs
