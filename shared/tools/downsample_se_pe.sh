#!/usr/bin/bash

## Usage: bash downsample_se_pe.sh 1000 10 a.fastq.gz a.downsampled.fastq.gz b.fastq.gz b.downsampled.fastq.gz

num_reads=$1
threads=$2
R1=$3
R1_out=$4
R2=$5
R2_out=$6

seed=12345

## check for paired or single end mode
paired=false
if [ -z "$R1" ] || [ -z "$R1_out" ]; then
	echo -e "\nPlease provide R1 input file and R1 output filename!\n"
	echo "usage: downsample_se_pe.sh 1000 10 a.fastq.gz a.downsampled.fastq.gz (single end)"
	echo "or"
	echo -e "(paired end) downsample_se_pe.sh 1000 10 a.fastq.gz a.downsampled.fastq.gz b.fastq.gz b.downsampled.fastq.gz\n"
	exit 1;
elif [ -n "$R2" ] && [ -n "$R2_out" ]; then
	paired=true
elif ( ( [ -z "$R2" ] && [ -n "$R2_out" ] ) || ( [ -n "$R2" ] && [ -z "$R2_out" ] ) ); then
	echo "Please provide for R2 either both input and output file or leave both empty!"
	exit 1;
fi

## check for pigz
tmp=$(which pigz)
if [ -z "$tmp" ]; then
	echo "Pigz needs to be installed!"
	exit 1;
fi

##############################################################################
## Usage: awk -v n=5 -v s=12345 -f downsample.awk <(seq 1 100)
## Reservoir downsampling in AWK!
## Fabian Kilpert, January 2017
read -d '' reservoir_sampling << 'EOF'
BEGIN {
  if (n=="") n = 1;

	if (s=="") srand()
	else srand(s)
}

# function generating numbers between 1 and n
function randint(n) { return int(1 + n * rand()) }

{
	# fill reservoir
	if (NR <= n) R[NR] = $0

	# replace reservoir elements with gradually decreasing probability
	else {
		j = randint(NR);
		if (j <= n) { R[j] = $0 };
	}
}

END {
	for (key in R) print R[key]
}
EOF
##############################################################################


## do the work
if [ "$paired" == true ]; then
	## paired end
	echo -e "\nDownsampling fastq files to $num_reads reads in paired-end mode..." > /dev/stderr;
	paste <(pigz -p2 -dc ${R1} | sed '/^$/d' | paste -d "\t" - - - - ) <(pigz -p2 -dc ${R2} | sed '/^$/d' | paste -d "\t" - - - - ) \
	| awk -v n=${num_reads} -v s=$seed "$reservoir_sampling" \
	| tee >(cut -f1,2,3,4 | sed 's/\t/\n/g' | pigz -p $(($threads/2>2?$threads/2:2)) -9 > ${R1_out}) >(cut -f5,6,7,8 | sed 's/\t/\n/g' | pigz -p $(($threads/2>2?$threads/2:2)) -9 > ${R2_out}) \
	| cut -f1,5 | tr " " "\t" | awk '{if ($1!=$3) print $0,"Read names in paired-end fastq files are out of sync!"; }'
else
	## single end
	echo -e "\nDownsampling fastq file to $num_reads reads in single-end mode..." > /dev/stderr
	pigz -p2 -dc ${R1} | sed '/^$/d' | paste -d "\t" - - - - | awk -v n=${num_reads} -v s=${seed} "$reservoir_sampling" | sed 's/\t/\n/g' | pigz -p $threads -9 > ${R1_out}
fi


################ END



#DIFF=$(diff <(pigz -dc ${2} | sed -n '1p;1~4p' | cut -d " " -f1) <(pigz -dc ${3} | sed -n '1p;1~4p' | cut -d " " -f1))
#if [ "$DIFF" != "" ] 
#then
#    echo "Error! FASTQ sequence identifiers do NOT match!"
#    exit 1
#else


