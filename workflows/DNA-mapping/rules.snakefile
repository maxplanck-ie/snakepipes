## FASTQ: downsample (on request) or create symlinks
if downsample:
    include: os.path.join(maindir, "shared", "rules", "FASTQ_downsample.snakefile")
else:
    include: os.path.join(maindir, "shared", "rules", "FASTQ_symlink.snakefile")


## FastQC
include: os.path.join(maindir, "shared", "rules", "FastQC.snakefile")


## TrimGalore
if trim:
    include: os.path.join(maindir, "shared", "rules", "TrimGalore.snakefile")


## Bowtie2_and_MarkDuplicates
##include: os.path.join(maindir, "shared", "rules", "Bowtie2.snakefile")
include: os.path.join(maindir, "shared", "rules", "Bowtie2_and_MarkDuplicates.snakefile")


## Picard CollectAlignmentSummaryMetrics and CollectInsertSizeMetrics
include: os.path.join(maindir, "shared", "rules", "Picard.snakefile")


## Qualimap bamqc
include: os.path.join(maindir, "shared", "rules", "Qualimap_bamqc.snakefile")
