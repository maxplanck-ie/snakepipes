mode=$1
outdir=$2
design=$3

dodesign=''
if [ "$design" != "" ] ; then
  dodesign="--DE $design"
fi

./RNA-seq -i ./sandbox/RNAseq_fastq/ -o $outdir mm10 -m $mode -j 10 --trim --fastqc --snakemake_options='--rerun-incomplete' $dodesign
