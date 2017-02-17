#!/bin/bash

dir_in=$1
dir_out=$2

dir_in=$(readlink -m $dir_in)

mkdir -p $dir_out/data

#for i in  $dir_in/*cout{b,c}.csv; do 
# out=$(echo $i | sed 's/.*\///'); 
# out2=$(echo $out | sed 's/\.csv$//');
#	echo $out 1>&2;
# cat $i | awk -v file=$out '{if (NR==1) {next;}; for (i=2;i<=NF;i++) COUNTS[i-1]+=$i;} \
#  END{for (i=1;i<=192;i++){if (i<=96) sum1+=COUNTS[i];else sum2+=COUNTS[i];} if (sum1>sum2) offset=1; else offset=97; \
#      for (i=offset;i<offset+96;i++) {OFS="\t";print file,i,COUNTS[i];sum+=COUNTS[i];}}' >$dir_out/data/$out.cellsum;
#done

if test -z "$(find $dir_in/ -maxdepth 1 -name '*.featureCounts_summary.txt')"; then 

  for i in  $dir_in/*.cout{b,c}.csv; do 
        type=$(echo $i | sed 's/.*\///' | sed 's/.*\.\(cout.\)\.csv$/\1/'); ## type = "coutc" or "coutb"
        sample=$(echo $i | sed 's/.*\///' | sed 's/\.cout.\.csv$//'); ## sample name without ending
	echo $out 1>&2;
	cat $i | awk -v sample=$sample type=$type '{
		if (NR==1) {cells = NF-1; next;}; 
		for (i=2;i<=NF;i++) COUNTS[i-1]+=$i;
	 }
 	END{
		#match(sample,"([^[:space:]\\.]+)\\.([^[:space:]\\.]+).csv",name)
		if (type~"coutc") print "sample\ncell_idx\tREADS_UNIQFEAT"; 
		else 	print "sample\ncell_idx\tUMI"; 	
		for (i=1;i<=cells;i++) {
			OFS="\t";print sample,i,COUNTS[i];
		}
	}' > $dir_out/data/$out.cellsum;
  done
else 
  
  for i in $dir_in/*.featureCounts_summary.txt; do
    out=$(cat $i | awk '{if (NR==2) print $1}');
    cat $i | sed -n -e '/sample.idx.READS/,/#LIB/{{/#LIB/d;p}}' > $dir_out/data/$out.cellsum;
  done


## with header!
paste <(echo -e "sample\tcell\tcell_reads"; cat $dir_out/data/*coutc*.cellsum | sort -k1,1 -k2,2n -V ) <(echo -e "sample2\tcell2\tcell_transcripts"; cat $dir_out/data/*coutb*.cellsum |  sort -k1,1 -k2,2n -V) > $dir_out/all_samples.cellsum_coutc_coutb.tsv

## summary per sample/lib
cat  $dir_out/all_samples.cellsum_coutc_coutb.tsv | awk '{if (NR<=2) last=$1; \
   if ($1==last) {sum1+=$3; sum2+=$6;} else {print last,sum1,sum2,sum2/sum1,sum1/96,sum2/96; sum1=0; sum2=0}; last=$1} \
   END{print last,sum1,sum2,sum2/sum1,sum1/96,sum2/96;}'

#> $dir/Results/QC_report.all_samples.tsv

#cat $dir/Results/QC_report.all_samples.tsv | column -t > $dir/Results/QC_report.all_samples.txt

