#!/bin/bash

## example call:
## scRNAseq_bam_featureCounts.sh test.sorted.bam genes.filtered.gtf celseq_barcodes.192.txt test.txt /package/subread-1.5.0-p1/bin/ tmp_fc 5 1>MySample.cout.csv 2>MySample.cout_summary.txt

bam=$1		## mapping for 96 cells
gtf=$2		## gene annotation
bc_file=$3		## celSeq cell barcode file
sample_name=$4	## sample name, used for featureCounts as ouput name, NOT a directory or path
lib_type=$5
tmp=$6	## will be created by this script, used as working dir for featureCounts due to -R issue
threads=$7	## for featureCOunts only


## gtf is expected in this format, we use only gene_id and gene_name
# 3	stdin	exon	108107280	108109316	.	-	.	gene_id "ENSMUSG00000000001.4"; transcript_id "ENSMUST00000000001.4"; exon_number "1"; exon_id "ENSMUST00000000001.4.1"; gene_name "Gnai3";
# 3	stdin	exon	108109403	108109612	.	-	.	gene_id "ENSMUSG00000000001.4"; transcript_id "ENSMUST00000000001.4"; exon_number "2"; exon_id "ENSMUST00000000001.4.2"; gene_name "Gnai3";

curr=$(pwd)

gtf_path=$(realpath $gtf)
bam_path=$(realpath $bam)
tmp_path=$(realpath $tmp)

mkdir -p $tmp_path 1>&2 
tmp_dir=$(mktemp -d --tmpdir=$tmp_path)

## call featureCounts
(>&2 echo "featureCounts -a $gtf_path -T $threads -s $lib_type -R CORE -F "GTF" --tmpDir ${tmp_dir} -o ${tmp_dir}/_tmp_$sample_name $bam_path 1>&2")
featureCounts -a $gtf_path -T $threads -s $lib_type -R CORE -F "GTF" --tmpDir ${tmp_dir} -o ${tmp_dir}/_tmp_$sample_name $bam_path 1>&2

## add gene_id (gtf col 10), gene_name (gtf col 18) to featureCounts output, last col is gene_name + chromosome
## <(cat $gtf_path | tr " " "\t" | tr -d "\";" | awk '{print $10,$18,$18"__chr"$1}') \
cat ${tmp_dir}/${sample_name}.bam.featureCounts | awk -v map_f=$gtf_path \
'BEGIN{
	while (getline < map_f) {
		match($0,"gene_id[[:space:]\";]+([^[:space:]\";]+)",gid)
		match($0,"gene_name[[:space:]\";]+([^[:space:]\";]+)",gna)
		MAP[gid[1]]=gna[1]"\t"gna[1]"__chr"$1; 
	} 
}
{OFS="\t";
if ($4 in MAP) print $0,MAP[$4]; else print $0,"NA","NA"; # MAP[$4] contains tab!
}' |
## Output after pipe is like:
## new featureCounts from 1.6
## J00182:58:HN77HBBXX:5:1115:11992:23681:SC:CGTTAC:41:UMI:TGGGGT:36:40:63 Unassigned_MultiMapping -1      NA	NA	NA
## J00182:58:HN77HBBXX:5:1213:15280:14150:SC:CGTTAC:41:UMI:TGGGGT:36:36:63 Unassigned_MultiMapping -1      NA	NA 	NA
## J00182:58:HN77HBBXX:5:2110:31507:27883:SC:AGTACC:41:UMI:ATCGGA:34:60:63 Assigned        1       ENSMUSG00000103201.1	Gm37180	Gm37180__chr1
## J00182:58:HN77HBBXX:6:2107:21724:18247:SC:AGTACC:41:UMI:ATCGGA:36:47:63 Assigned        1       ENSMUSG00000103201.1	Gm37180	Gm37180__chr1
##
## old featurecounts <=1.5.2
## SN7001180:281:C99CMACXX:2:1203:1365:44875:SC:ACTCGA:37:UMI:ATTCCT:35:36:38	Unassigned_NoFeatures   *       *	NA	NA
## SN7001180:281:C99CMACXX:2:2212:8564:84823:SC:ACGTGA:37:UMI:CGCCAG:35:35:38	Assigned	ENSMUSG00000103377.1	*	Gm37180	Gm37180__chr1
## SN7001180:281:C99CMACXX:2:1211:18877:23349:SC:GACAAC:37:UMI:TGTCCG:35:27:38	Unassigned_Ambiguity	*	Number_Of_Overlapped_Genes=2	NA	NA
##
## now we can count by getting the cell barcode and UMI from readname
## put all in big matrix in awk and write to stdout to caputure this later
## summary stats are printed to stderr
##
awk -v map_f=$bc_file -v sample=$sample_name ' \
BEGIN{
	while(getline<map_f) {                      ## read in cell barcodes
		CELL[$2]=$1; num_cells+=1;
	}
	nocell_unmap=0;
	nocell_map=0;
}
{
	if ($1 in READS_SEEN) next;
	
	pos=match($1,":SC:");                       ## get barcode startpos (":SC:"") from readname
	split(substr($1,pos+1),BC,":");             ## split on ":" to separate all info and stor in array "BC"

	if (BC[2] in CELL) {					## valid cell barcode
		if (BC[5]~"N") cell_noumi[CELL[BC[2]]] += 1;    ## UMIs with N are also not considered later on
		else if ($2~"Assigned" && $3 != "-1") {         ## if unique feature then count it
			ALL[$6][BC[5]][CELL[BC[2]]] += 1;         ##
			cell_uniqfeat[CELL[BC[2]]] += 1; }        ## only stats
		else if ($2~"NoFeatures") cell_nofeat[CELL[BC[2]]] += 1;
		else if ($2~"MultiMapping") cell_multimap[CELL[BC[2]]] += 1;			
		else if ($2~"Unassigned_Ambiguity") cell_multifeat[CELL[BC[2]]] +=1;
		else if ($2~"Unassigned_Unmapped") cell_unmap[CELL[BC[2]]] +=1 ;
	} else if ($2~"Unassigned_Unmapped") nocell_unmap+=1; 
	else nocell_map+=1;

	if ($2!~"Unassigned_Unmapped") READS_SEEN[$1];	## only for unmapped reads it is safe to ignore this check 
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

	print "sample\tcell_idx\tREADS_UNMAP\tREADS_NOFEAT\tREADS_NOUMI\tREADS_MULTIMAP\tREADS_MULTIFEAT\tREADS_UNIQFEAT\tUMI" > "/dev/stderr";
	for (j=1;j<=num_cells;j++) {
		out = sample"\t"j;
		if ( j in cell_unmap) out = out"\t"cell_unmap[j]; else out = out"\t0";
		if ( j in cell_nofeat) out = out"\t"cell_nofeat[j]; else out = out"\t0";
		if ( j in cell_noumi) out = out"\t"cell_noumi[j]; else out = out"\t0";
		if ( j in cell_multimap) out = out"\t"cell_multimap[j]; else out = out"\t0";
		if ( j in cell_multifeat) out = out"\t"cell_multifeat[j]; else out = out"\t0";
		if ( j in cell_uniqfeat) out = out"\t"cell_uniqfeat[j]; else out = out"\t0";
		if ( j in cell_UMI) out = out"\t"cell_UMI[j]; else out = out"\t0";
		print out > "/dev/stderr";

		ALLcell_unmap += cell_unmap[j];
		ALLcell_nofeat += cell_nofeat[j];
        ALLcell_noumi += cell_noumi[j];
		ALLcell_multimap += cell_multimap[j];
		ALLcell_multifeat += cell_multifeat[j];
		ALLcell_uniqfeat += cell_uniqfeat[j];
		ALLcell_UMI += cell_UMI[j];
	}

	sum_reads = ALLcell_uniqfeat + ALLcell_nofeat + ALLcell_noumi + ALLcell_multifeat + ALLcell_multimap + nocell_map;
	sum_reads_all = sum_reads + nocell_unmap + ALLcell_unmap;
	sum = "#LIBREADS_UMI\t"ALLcell_UMI"\t"(ALLcell_UMI/sum_reads*100)"\n";
	sum = sum"#LIBREADS_UNIQFEAT\t"ALLcell_uniqfeat"\t"(ALLcell_uniqfeat/sum_reads*100)"\n";
	sum = sum"#LIBREADS_MULTIMAP\t"ALLcell_multimap"\t"(ALLcell_multimap/sum_reads*100)"\n";
	sum = sum"#LIBREADS_MULTIFEAT\t"ALLcell_multifeat"\t"(ALLcell_multifeat/sum_reads*100)"\n";
	sum = sum"#LIBREADS_NOUMI\t"ALLcell_noumi"\t"(ALLcell_noumi/sum_reads*100)"\n";	
	sum = sum"#LIBREADS_NOFEAT\t"ALLcell_nofeat"\t"(ALLcell_nofeat/sum_reads*100)"\n";
	sum = sum"#LIBREADS_NOCELL\t"nocell_map"\t"(nocell_map/sum_reads*100)"\n";
	sum = sum"#LIBREADS_MAPTOTAL\t"sum_reads"\t100.0\n";
	sum = sum"#LIBREADS_UNMAPCELL\t"ALLcell_unmap"\t"(ALLcell_unmap/sum_reads_all*100)"\n";
	sum = sum"#LIBREADS_UNMAP\t"nocell_unmap"\t"(nocell_unmap/sum_reads_all*100)"\n";
	sum = sum"#LIBREADS_TOTAL\t"sum_reads_all"\t100.0";
	print sum > "/dev/stderr";                  ## prints stats to stderr
}' 

#2> >(tee >(grep "^#" | tr -d "#" > test_sum.txt) >(grep -v "^#" > test_cell.txt))

