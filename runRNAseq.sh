mode=$1
outdir=$2
design=$3

FASTQ=/data/manke/group/rauer/pospisilikgroup/BimodalMice/data/Project_A802/RNAseq_fastq

dodesign=''
if [ "$design" != "" ] ; then
  dodesign="--DE $design"
fi



./RNA-seq -i $FASTQ -o $outdir mm10 -m $mode -j 10 --trim --fastqc --snakemake_options='--rerun-incomplete' $dodesign
