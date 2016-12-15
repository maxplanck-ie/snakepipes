#!/bin/bash

bam=$1
gtf=$2
bc_file=$3
out=$4
fc_path=$5
tmp=$6

curr=$(pwd)

gtf_path=$(realpath $gtf)
bam_path=$(realpath $bam)

mkdir -p $tmp
cd $tmp

${fc_path}featureCounts -a $gtf_path -s 1 -R -d 25 -F "GTF" -o $out $bam_path 

cat *.bam.featureCounts | awk -v map_f=<(cat $gtf_path | tr " " "\t" | tr -d "\";" | awk '{print $10,$18,$18"__chr"$1}') \
'BEGIN{while (getline < map_f) { MAP[$1]=$2"\t"$3; } }
{OFS="\t";
if ($3 in MAP) print $0,MAP[$3]; else print $0,"NA","NA";
}' |
#END {print "noFEAT",nofeat; print "Ass",ass_feat; print "Amb",feat_amb;}' 
## SN7001180:281:C99CMACXX:2:1203:1365:44875:SC:ACTCGA:37:UMI:ATTCCT:35:36:38      Unassigned_NoFeatures   *       *
## SN7001180:281:C99CMACXX:2:1306:16709:23741:SC:AGCTAG:37:UMI:TGGAGA:35:22:38     Assigned        ENSMUSG00000025907.14   *
## SN7001180:281:C99CMACXX:2:2302:19826:45553:SC:AGCTAG:37:UMI:TGGAGA:35:34:38     Assigned        ENSMUSG00000025907.14   *
## SN7001180:281:C99CMACXX:2:1211:18877:23349:SC:GACAAC:37:UMI:TGTCCG:35:27:38	Unassigned_Ambiguity	*	Number_Of_Overlapped_Genes=2
## SN7001180:281:C99CMACXX:2:2202:2909:75663:SC:GACAAC:37:UMI:TGTCCG:35:30:38	Unassigned_Ambiguity	*	Number_Of_Overlapped_Genes=2
awk -v map_f=$bc_file ' \
BEGIN{
	while(getline<map_f) {                      ## read in cell barcodes
		CELL[$2]=$1; num_cells+=1;
	}
}
{
	pos=match($1,":SC:");                       ## get startpos of precious barcode info from readname
	split(substr($1,pos+1),BC,":");             ## split on ":" to separate all info and stor in array "BC"
	num=split($3,GENES,",");                    ## get number of features a read overlaps by splitting on ","
	if ( $3!="*" && num==1 && BC[2] in CELL) {  ## if unique feature and a "valid" cell barcode then count it
		ALL[$6][BC[5]][CELL[BC[2]]] += 1;     ## $3 is single GENEID if num==1
		feat_uniq+=1;                         ## only stats
	}
	if (BC[2] in CELL) {                        ## only stats
		if ($2~"NoFeatures") cell_nofeat+=1;
			else 
			 if ($2~"Unassigned_Ambiguity") cell_multi+=1; else if ($2~"Assigned") feat_uniq1+=1;    
	} else nocell+=1;
}
END{
	printf "GENEID\tRBAR";                      ## mimic Dominics output format
	for (n=1;n<=num_cells;n++)                  ## header line  
		printf "\t"n;
	printf "\n";
	for (i in ALL) {                            ## iterate over all genes
		for (k in ALL[i]) {                   ## and all UMIs
			printf i"\t"k;
			for (j=1;j<=num_cells;j++) {    ## and all cells
				if (j in ALL[i][k])
					printf "\t"ALL[i][k][j];
				else printf "\t0";
			}
			printf "\n";
		}
	}
	sum_reads = feat_uniq1 + cell_nofeat + cell_multi + nocell;
	sum = "FEATURE_UNIQUE\t"feat_uniq1"\t"(feat_uniq1/sum_reads*100)"\n";
	sum = sum"FEATURE_MULTI\t"cell_multi"\t"(cell_multi/sum_reads*100)"\n";
	sum = sum"CELL_NOFEATURE\t"cell_nofeat"\t"(cell_nofeat/sum_reads*100)"\n";
	sum = sum"NOCELL\t"nocell"\t"(nocell/sum_reads*100)"\n";
	sum = sum"NUM_READS\t"sum_reads"\t100.0";
	print sum > "/dev/stderr";                  ## prints stats to stderr
}' 

cd $curr
