#! /usr/bin/env bash

GTF=$1
OUTDIR=$2

T2G=$OUTDIR/$(basename ${GTF%%".gtf"}".t2g")
# GTF=/data/repository/organisms/GRCm38_ensembl/gencode/m9/genes.gtf
cat $GTF | grep -v '^##' | head -n2000 | awk '{ if ( $3 == "transcript" ) {print ($0) }}' | cut -f 9 | awk -F ";" 'BEGIN {OFS=";"}{  print ($2,$1,$5) } ' > tmp

for key in gene_id transcript_id gene_name "\"" " ";do  sed -i -e "s/$key//g" tmp ;done;

sed -i -e "s/;/\t/g" tmp

echo -e "gene_id\ttranscript_id\tgene_name" > $T2G
cat tmp >> $T2G
