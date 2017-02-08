#!/bin/bash

## example call:
## scRNAseq_bam_featureCounts.sh test.sorted.bam genes.filtered.gtf celseq_barcodes.192.txt test.txt /package/subread-1.5.0-p1/bin/ tmp_fc 5 1>MySample.cout.csv 2>MySample.cout_summary.txt

bam=$1		## mapping for 96 cells
gtf=$2		## gene annotation
bc_file=$3	## celSeq cell barcode file
out=$4		## just used for featureCounts as ouput name, NOT a directory or path
fc_path=$5	## path to fc like "/package/subread-1.5.0-p1/bin/"
tmp=$6		## will be created by this script, used as working dir for featureCounts due to -R issue
threads=$7	## for featureCOunts only


## gtf is expected in this format, we use only gene_id and gene_name
# 3	stdin	exon	108107280	108109316	.	-	.	gene_id "ENSMUSG00000000001.4"; transcript_id "ENSMUST00000000001.4"; exon_number "1"; exon_id "ENSMUST00000000001.4.1"; gene_name "Gnai3";
# 3	stdin	exon	108109403	108109612	.	-	.	gene_id "ENSMUSG00000000001.4"; transcript_id "ENSMUST00000000001.4"; exon_number "2"; exon_id "ENSMUST00000000001.4.2"; gene_name "Gnai3";

curr=$(pwd)

gtf_path=$(realpath $gtf)
bam_path=$(realpath $bam)
tmp_path=$(realpath $tmp)

## current version of featureCounts under /package/subread-1.5.0-p1/ writes out -R file to currDir instead to path provided with -o
## this is fixed in more recent version of subread! We have to install it! :-)
mkdir -p $tmp_path
cd $tmp_path
rm *.bam.featureCounts		## I'm too lazy the get the full correct name later on, so make sure we have only the file we want

## call featureCounts
${fc_path}featureCounts -a $gtf_path -T $threads -s 1 -R -d 25 -F "GTF" -o _tmp_$out $bam_path 

## add gene_id (gtf col 10), gene_name (gtf col 18) to featureCounts output, last col is gene_name + chromosome
## <(cat $gtf_path | tr " " "\t" | tr -d "\";" | awk '{print $10,$18,$18"__chr"$1}') \
cat *.bam.featureCounts | awk -v map_f=$gtf_path \
'BEGIN{
	while (getline < map_f) {
		match($0,"gene_id[[:space:]\";]+([^[:space:]\";]+)",gid)
		match($0,"gene_name[[:space:]\";]+([^[:space:]\";]+)",gna)
		MAP[gid[1]]=gna[1]"\t"gna[1]"__chr"$1; 
	} 
}
{OFS="\t";
if ($3 in MAP) print $0,MAP[$3]; else print $0,"NA","NA";
}' |
## Output after pipe is like:
##
## SN7001180:281:C99CMACXX:2:1203:1365:44875:SC:ACTCGA:37:UMI:ATTCCT:35:36:38	Unassigned_NoFeatures   *       *	NA	NA
## SN7001180:281:C99CMACXX:2:2212:8564:84823:SC:ACGTGA:37:UMI:CGCCAG:35:35:38	Assigned	ENSMUSG00000103377.1	*	Gm37180	Gm37180__chr1
## SN7001180:281:C99CMACXX:2:1211:18877:23349:SC:GACAAC:37:UMI:TGTCCG:35:27:38	Unassigned_Ambiguity	*	Number_Of_Overlapped_Genes=2	NA	NA
##
## now we can count by getting the cell barcode and UMI from readname
## put all in big matrix in awk and write to stdout to caputure this later
## summary stats are printed to stderr
##
awk -v map_f=$bc_file ' \
BEGIN{
	while(getline<map_f) {                      ## read in cell barcodes
		CELL[$2]=$1; num_cells+=1;
	}
}
{
	pos=match($1,":SC:");                       ## get barcode startpos (":SC:"") from readname
	split(substr($1,pos+1),BC,":");             ## split on ":" to separate all info and stor in array "BC"
#	if ( $3!="*" && BC[2] in CELL) {  		## if unique feature and a "valid" cell barcode then count it
#		ALL[$6][BC[5]][CELL[BC[2]]] += 1;   ## $3 is single GENEID if num==1
#		feat_uniq+=1;                       ## only stats
#	}
	if (BC[2] in CELL) {
		if ($2~"Assigned" && $3 != "*") {
			ALL[$6][BC[5]][CELL[BC[2]]] += 1;
			cell_uniqfeat[CELL[BC[2]]] += 1; }    
		else if ($2~"NoFeatures") cell_nofeat[CELL[BC[2]]] += 1;
		else if ($2~"Multimapping") cell_multimap[CELL[BC[2]]] += 1;
		else if ($2~"Unassigned_Ambiguity") cell_multifeat[CELL[BC[2]]] +=1 ;
	} else nocell+=1;
}
END{
	printf "GENEID\tRBAR";                        ## mimic Dominics output format
	for (n=1;n<=num_cells;n++)                    ## header line  
		printf "\t"n;
	printf "\n";
	for (i in ALL) {                              ## iterate over all genes
		for (k in ALL[i]) {                   ## and all UMIs
			printf i"\t"k;
			for (j=1;j<=num_cells;j++) {  ## and all cells
				if (j in ALL[i][k]){
					printf "\t"ALL[i][k][j];
					cell_UMI[j] += 1;
				}
				else printf "\t0";
			}
			printf "\n";
		}
	}

	print "#idx\tNOFEAT\tMULTIMAP\tMULTIFEAT\tUNIQFEAT\tUMI" > "/dev/stderr";
	for (j=1;j<=num_cells;j++) {
		out = ""j;
		if ( j in cell_nofeat) out = out"\t"cell_nofeat[j]; else out = out"\t0";
		if ( j in cell_multimap) out = out"\t"cell_multimap[j]; else out = out"\t0";
		if ( j in cell_multifeat) out = out"\t"cell_multifeat[j]; else out = out"\t0";
		if ( j in cell_uniqfeat) out = out"\t"cell_uniqfeat[j]; else out = out"\t0";
		if ( j in cell_UMI) out = out"\t"cell_UMI[j]; else out = out"\t0";
		print out > "/dev/stderr";

		ALLcell_nofeat += cell_nofeat[j];
		ALLcell_multimap += cell_multimap[j];
		ALLcell_multifeat += cell_multifeat[j];
		ALLcell_uniqfeat += cell_uniqfeat[j];
		ALLcell_UMI += cell_UMI[j];
	}

	sum_reads = ALLcell_uniqfeat + ALLcell_nofeat + ALLcell_multifeat + ALLcell_multimap + nocell;
	sum = "#FEATURE_UNIQUE\t"ALLcell_uniqfeat"\t"(ALLcell_uniqfeat/sum_reads*100)"\n";
	sum = sum"#UMI\t"ALLcell_UMI"\t"(ALLcell_UMI/sum_reads*100)"\n";
	sum = sum"#READS_MULTIMAP\t"ALLcell_multimap"\t"(ALLcell_multimap/sum_reads*100)"\n";
	sum = sum"#FEATURE_MULTI\t"ALLcell_multifeat"\t"(ALLcell_multifeat/sum_reads*100)"\n";
	sum = sum"#CELL_NOFEATURE\t"ALLcell_nofeat"\t"(ALLcell_nofeat/sum_reads*100)"\n";
	sum = sum"#NOCELL\t"nocell"\t"(nocell/sum_reads*100)"\n";
	sum = sum"#NUM_READS\t"sum_reads"\t100.0";
	print sum > "/dev/stderr";                  ## prints stats to stderr
}'

