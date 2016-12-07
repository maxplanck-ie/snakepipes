#!/bin/bash

dir=$1

dir=$(realpath $dir)

mkdir -p $dir/Results/QC

for i in  $dir/Counts/*cout{b,c}.csv; do 
 out=$(echo $i | sed 's/.*\///'); 
	echo $out;
 cat $i | awk -v file=$out '{if (NR==1) {next;}; for (i=2;i<=NF;i++) COUNTS[i-1]+=$i;} \
  END{for (i=1;i<=192;i++){if (i<=96) sum1+=COUNTS[i];else sum2+=COUNTS[i];} if (sum1>sum2) offset=1; else offset=97; \
      for (i=offset;i<offset+96;i++) {OFS="\t";print file,i,COUNTS[i];sum+=COUNTS[i];}}' >$dir/Results/QC/$out.cellsum;
done

paste <(cat $dir/Results/QC/*coutc*.cellsum) <(cat $dir/Results/QC/*coutb*.cellsum) > $dir/Results/all_samples.cellsum_coutc_countb.tsv

cat  $dir/Results/all_samples.cellsum_coutc_countb.tsv | awk '{if (NR==1) last=$1; \
   if ($1==last) {sum1+=$3; sum2+=$6;} else {print last,sum1,sum2,sum2/sum1,sum1/96,sum2/96; sum1=0; sum2=0}; last=$1} \
   END{print last,sum1,sum2,sum2/sum1,sum1/96,sum2/96;}' > $dir/Results/QC_report.all_samples.tsv

cat $dir/Results/QC_report.all_samples.tsv | column -t > $dir/Results/QC_report.all_samples.txt

