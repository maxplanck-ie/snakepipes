#!/bin/bash

bam=$1
gtf=$2
out=$3
fc_path=$4
tmp=$5

curr=$(pwd)
cd $tmp

${fc_path}featureCounts -a $gtf -s 1 -R -d 25 -F "GTF" -o $out $bam


