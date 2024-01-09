#!/usr/bin/env bash

set -ue

if [ $# -lt 5 ]
then
	echo "
	SEACR: Sparse Enrichment Analysis for CUT&RUN
	
	Usage: bash SEACR_1.3.sh <experimental bedgraph>.bg [<control bedgraph>.bg | <FDR threshold>] ["norm" | "non"] ["relaxed" | "stringent"] output prefix
	
	Description of input fields:
	
	Field 1: Target data bedgraph file in UCSC bedgraph format (https://genome.ucsc.edu/goldenpath/help/bedgraph.html) that omits regions containing 0 signal.
	
	Field 2: Control (IgG) data bedgraph file to generate an empirical threshold for peak calling. Alternatively, a numeric threshold n between 0 and 1 returns the top n fraction of peaks based on total signal within peaks.
	
	Field 3: “norm” denotes normalization of control to target data, “non” skips this behavior. "norm" is recommended unless experimental and control data are already rigorously normalized to each other (e.g. via spike-in).
		
	Field 4: “relaxed” uses a total signal threshold between the knee and peak of the total signal curve, and corresponds to the “relaxed” mode described in the text, whereas “stringent” uses the peak of the curve, and corresponds to “stringent” mode.
	
	Field 5: Output prefix
	
	Output file:
	<output prefix>.auc.threshold.merge.bed (Bed file of enriched regions)
	
	Output data structure: 
	
	<chr>	<start>	<end>	<AUC>	<max signal>	<max signal region>
	
	Description of output fields:
	Field 1: Chromosome
	
	Field 2: Start coordinate
	
	Field 3: End coordinate
	
	Field 4: Total signal contained within denoted coordinates
	
	Field 5: Maximum bedgraph signal attained at any base pair within denoted coordinates
	
	Field 6: Region representing the farthest upstream and farthest downstream bases within the denoted coordinates that are represented by the maximum bedgraph signal
	
	Examples:
	bash SEACR_1.3.sh target.bedgraph IgG.bedgraph norm stringent output
	Calls enriched regions in target data using normalized IgG control track with stringent threshold
	
	bash SEACR_1.3.sh target.bedgraph IgG.bedgraph non relaxed output
	Calls enriched regions in target data using non-normalized IgG control track with relaxed threshold
	bash SEACR_1.3.sh target.bedgraph 0.01 non stringent output
	Calls enriched regions in target data by selecting the top 1% of regions by area under the curve (AUC)
	"
	exit 1
fi

password=`head /dev/urandom | LC_CTYPE=C tr -dc A-Za-z0-9 | head -c 13; echo ''`
password2=`head /dev/urandom | LC_CTYPE=C tr -dc A-Za-z0-9 | head -c 13; echo ''`

exp=`basename $1`

if [[ $2 =~ ^[0-9]?+([.][0-9]+)?$ ]] || [[ $2 =~ ^[0-9]([.][0-9]+) ]] || [[ $2 =~ ^([.][0-9]+) ]]
then
	echo "Calling enriched regions without control file"
elif [[ -f $2 ]]
then
	echo "Calling enriched regions with control file"
	ctrl=`basename $2`
else
	echo "$2 is not a number or a file"
	exit 1
fi

norm=`echo $3`

if [[ $norm == "norm" ]]
then
	echo "Normalizing control to experimental bedgraph"
elif [[ $norm == "non" ]]
	then
	echo "Proceeding without normalization of control to experimental bedgraph"
else
	echo "Must specify \"norm\" for normalized or \"non\" for non-normalized data processing in third input"
	exit 1
fi

height=`echo $4`

if [[ $height == "relaxed" ]]
then
	echo "Using relaxed threshold"
elif [[ $height == "stringent" ]]
	then
	echo "Using stringent threshold"
else
	echo "Must specify \"stringent\" or \"relaxed\" in fourth input"
	exit 1
fi

echo "Creating experimental AUC file: $(date)"

awk 'BEGIN{s=1}; {if(s==1){s++}else if(s==2){if($4 > 0){chr=$1; start=$2; stop=$3; max=$4; coord=$1":"$2"-"$3; auc=$4*($3-$2); num=1; s++}}else{if($4 > 0){if(chr==$1 && $2==stop){num++; stop=$3; auc=auc+($4*($3-$2)); if ($4 > max){max=$4; coord=$1":"$2"-"$3
}else if($4 == max){split(coord,t,"-"); coord=t[1]"-"$3}}else{print chr"\t"start"\t"stop"\t"auc"\t"max"\t"coord"\t"num; chr=$1; start=$2; stop=$3; max=$4; coord=$1":"$2"-"$3; auc=$4*($3-$2); num=1}}}}' $1 > $password.auc.bed
cut -f 4,7 $password.auc.bed > $password.auc

if [[ -f $2 ]]
then
  echo "Creating control AUC file: $(date)"

  awk 'BEGIN{s=1}; {if(s==1){s++}else if(s==2){if($4 > 0){chr=$1; start=$2; stop=$3; max=$4; coord=$1":"$2"-"$3; auc=$4*($3-$2); num=1; s++}}else{if($4 > 0){if(chr==$1 && $2==stop){num++; stop=$3; auc=auc+($4*($3-$2)); if ($4 > max){max=$4; coord=$1":"$2"-"
$3}else if($4 == max){split(coord,t,"-"); coord=t[1]"-"$3}}else{print chr"\t"start"\t"stop"\t"auc"\t"max"\t"coord"\t"num; chr=$1; start=$2; stop=$3; max=$4; coord=$1":"$2"-"$3; auc=$4*($3-$2); num=1}}}}' $2 > $password2.auc.bed
  cut -f 4,7 $password2.auc.bed > $password2.auc
fi

# module load R  ## For use on cluster

echo "Calculating optimal AUC threshold: $(date)"

path=`dirname $0`
if [[ -f $2 ]] && [[ $norm == "norm" ]]
then
	echo "Calculating threshold using normalized control: $(date)"
	Rscript $path/SEACR_1.3.R --exp=$password.auc --ctrl=$password2.auc --norm=yes --output=$password
elif [[ -f $2 ]]
then
	echo "Calculating threshold using non-normalized control: $(date)"
	Rscript $path/SEACR_1.3.R --exp=$password.auc --ctrl=$password2.auc --norm=no --output=$password
else
	echo "Using user-provided threshold: $(date)"
	Rscript $path/SEACR_1.3.R --exp=$password.auc --ctrl=$2 --norm=no --output=$password
fi
	
fdr=`cat $password.fdr.txt | sed -n '1p'`			## Added 5/15/19 for SEACR_1.1
fdr2=`cat $password.fdr.txt | sed -n '2p'`			## Added 5/15/19 for SEACR_1.1

#thresh=`cat $exp.threshold.txt`
thresh=`cat $password.threshold.txt | sed -n '1p'`
thresh2=`cat $password.threshold.txt | sed -n '2p'`
thresh3=`cat $password.threshold.txt | sed -n '3p'`

echo "Creating thresholded feature file: $(date)"

if [[ $height == "relaxed" ]]
then
  echo "Empirical false discovery rate = $fdr2"
  awk -v value=$thresh2 -v value2=$thresh3 '$4 > value && $7 > value2 {print $0}' $password.auc.bed | cut -f 1,2,3,4,5,6 > $password.auc.threshold.bed
else
  echo "Empirical false discovery rate = $fdr"
  awk -v value=$thresh -v value2=$thresh3 '$4 > value && $7 > value2 {print $0}' $password.auc.bed | cut -f 1,2,3,4,5,6 > $password.auc.threshold.bed
fi

if [[ -f $2 ]]
then
	if [[ $norm == "norm" ]] #If normalizing, multiply control bedgraph by normalization constant
	then
		constant=`cat $password.norm.txt | sed -n '1p'`
		awk -v mult=$constant 'BEGIN{OFS="\t"}; {$4=$4*mult; print $0}' $password2.auc.bed | cut -f 1,2,3,4,5,6 > $password2.auc2.bed
		mv $password2.auc2.bed $password2.auc.bed
	fi
	awk -v value=$thresh '$4 > value {print $0}' $password2.auc.bed > $password2.auc.threshold.bed
fi

echo "Merging nearby features and eliminating control-enriched features: $(date)"

# module load bedtools ## For use on cluster
mean=`awk '{s+=$3-$2; t++}END{print s/(t*10)}' $password.auc.threshold.bed`

if [[ -f $2 ]]
then
	awk -v value=$mean 'BEGIN{s=1}; {if(s==1){chr=$1; start=$2; stop=$3; auc=$4; max=$5; coord=$6; s++}else{if(chr==$1 && $2 < stop+value){stop=$3; auc=auc+$4; if($5 > max){max=$5; coord=$6}else if($5==max){split(coord,t,"-"); split($6,u,"-"); coord=t[1]"-"u[2]}}else{print chr"\t"start"\t"stop"\t"auc"\t"max"\t"coord; chr=$1; start=$2; stop=$3; auc=$4; max=$5; coord=$6}}}' $password.auc.threshold.bed | bedtools intersect -wa -v -a - -b $password2.auc.threshold.bed > $5.auc.threshold.merge.bed  
else
	awk -v value=$mean 'BEGIN{s=1}; {if(s==1){chr=$1; start=$2; stop=$3; auc=$4; max=$5; coord=$6; s++}else{if(chr==$1 && $2 < stop+value){stop=$3; auc=auc+$4; if($5 > max){max=$5; coord=$6}else if($5==max){split(coord,t,"-"); split($6,u,"-"); coord=t[1]"-"u[2]}}else{print chr"\t"start"\t"stop"\t"auc"\t"max"\t"coord; chr=$1; start=$2; stop=$3; auc=$4; max=$5; coord=$6}}}' $password.auc.threshold.bed > $5.auc.threshold.merge.bed
fi

if [[ $height == "relaxed" ]]
then
  cat $5.auc.threshold.merge.bed > $5.relaxed.bed
else
  cat $5.auc.threshold.merge.bed > $5.stringent.bed
fi

echo "Removing temporary files: $(date)"

rm $password.auc.bed
rm $password.auc
rm $password.threshold.txt
rm $password.auc.threshold.bed
rm $password.fdr.txt  ## Added 5/15/19 for SEACR_1.1
rm $5.auc.threshold.merge.bed
if [[ -f $2 ]]
then
	rm $password2.auc.bed
	rm $password2.auc
	rm $password2.auc.threshold.bed
fi
if [[ $norm == "norm" ]]
then
	rm -f $password.norm.txt
fi
echo "Done: $(date)"
