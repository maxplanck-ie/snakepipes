#!/bin/bash

## example call:
## scRNAse_bam_genomic_feature_count.sh B6_11wks_1.bam genes.filtered.bed celseq_barcodes.192.txt /package/bedtools2-2.25.0/bin/ /package/samtools/bin/

bam=$1
bed=$2
bc_file=$3
bedtools=$4
samtools=$5

## -a is bam with genomic mapped reads with barcodes in read name
## -b is just bed file of features, but we want it sorted exactly as bam @SQ entries and also skip features on not-in-bam chr
## same sorting of -a -b files speeds up intersectBed as it can be used with "-sorted" option!
## hence we use 'sortBed' with faidx from bam file and use awk to skip bed entries on chr. not in bam
## -g is again faidx created from bam for 'sorted' option
${bedtools}intersectBed \
-a $bam \
-b <(${bedtools}sortBed -faidx <(${samtools}samtools view -H $bam | grep @SQ | sed 's/SN:\|LN://g' | awk '{OFS="\t";print $2,$3}') \
                        -i <(cat $bed | awk -v faidx_f=<(${samtools}samtools view -H $bam | grep @SQ | sed 's/SN:\|LN://g') \
                                         'BEGIN{while (getline < faidx_f) IDX[$2]=$2;}{OFS="\t";if ($1 in IDX) print $0;}' \
                                      | cut -f1-12 \
                            ) \
    ) \
-g <(${samtools}samtools view -H $bam | grep @SQ | sed 's/SN:\|LN://g' | awk '{OFS="\t";print $2,$3}') \
-sorted -split -bed -wao -s -nonamecheck |
## output at this point looks like: cols 1-12 from bam; 13-24 from bed, 25 #nt_overlap
## 1	3361041	3361087	SN7001180:281:C99CMACXX:2:2201:18880:15167:SC:ATCACG:37:UMI:TGGATC:35:30:38	60	-	3361041	3361087	0,0,0	1	46,	0,	.	-1	-1	.	-1	.	..	.	.	.	.	0
## 1	3667588	3667637	SN7001180:281:C99CMACXX:2:2201:15003:5111:SC:CTGTTG:37:UMI:AACCGT:35:25:38	60	+	3667588	3667637	0,0,0	1	49,	0,	.	-1	-1	.	-1	.	..	.	.	.	.	0
## 1	4491569	4491619	SN7001180:281:C99CMACXX:2:2201:15070:4540:SC:GTGACA:37:UMI:ATCTTA:35:30:38	60	-	4491569	4491619	0,0,0	1	50,	0,	1	4490930	4496413	ENSMUST00000027035.9	0	-	4491715	4493406	0	5	1738,367,92,807,123,	0,2169,2841,4205,5360,	50
## 1	4491569	4491619	SN7001180:281:C99CMACXX:2:2201:15070:4540:SC:GTGACA:37:UMI:ATCTTA:35:30:38	60	-	4491569	4491619	0,0,0	1	50,	0,	1	4491249	4496757	ENSMUST00000195555.1	0	-	4491715	4492591	0	3	1414,96,467,	0,2518,5041,	50
## 1	4491569	4491619	SN7001180:281:C99CMACXX:2:2201:15070:4540:SC:GTGACA:37:UMI:ATCTTA:35:30:38	60	-	4491569	4491619	0,0,0	1	50,	0,	1	4491389	4497354	ENSMUST00000192650.
##
## we use awk to add gene_name and gene_name_chr, we get this info from bed file columns: 13(transcript_id),15(gene_id),16(gene_name)
## new col 26: gene-id mapped from pipe output col 16 (transcript-id)
## new col 27: gene_name_chr
awk -v map_f=$bed 'BEGIN{while (getline < map_f) {MAP[$13]=$15;MAP2[$13]=$16} } 
{OFS="\t";if ($13!=".") print $0,MAP[$16],MAP2[$16]"__chr"$1; else print $0,"NA","NA"; } ' |
## 
## output is like:
##
## 1       3647281 3647322 SN7001180:281:C99CMACXX:2:1310:5055:22323:SC:AACACA:37:UMI:GGCTTC:35:14:38      60      +       3647281 3647322 0,0,0   1       41,     0,      .       -1      -1      .       -1      . .       .       .       .       .       .       0       NA      NA
## 1       3647310 3647357 SN7001180:281:C99CMACXX:2:2201:10782:57572:SC:AACCTC:34:UMI:AGACCG:35:38:38     60      -       3647310 3647357 0,0,0   1       47,     0,      1       3647308 3658904 ENSMUST00000192692.1    0       -       3658904 3658904 0       2       3201,58,        0,11538,        47      ENSMUSG00000102331.1    Gm19938__chr1
## 1       3647314 3647362 SN7001180:281:C99CMACXX:2:1214:18362:36904:SC:TGTCGA:37:UMI:GTTGTT:35:35:38     60      -       3647314 3647362 0,0,0   1       48,     0,      1       3647308 3658904 ENSMUST00000192692.1    0       -       3658904 3658904 0       2       3201,58,        0,11538,        48      ENSMUSG00000102331.1    Gm19938__chr1
## 1
##
## now we collect all overalpping features by read name, we keep distinct gene-ids, gene-names per read as comma-sep. list
## col 2 are all overlapping gene-ids for that read, col3 the corresponding gene_name_chr and col4 all transcript-ids, 
## last sed is to get rid of "," at the end of the collected ID lists
##
awk '{
	if ($4!=last) { 
		if (NR!=1 && no_ol==0) { 
			printf last"\t"; 
			for (i in GID) printf i","; 
			printf"\t";
			for (i in GNA) printf i",";
			printf "\t"; 
			for (i in TID) printf i","; 
			printf"\n"; 
			delete GID; delete TID; delete GNA;
		} else if (no_ol==1)
			print last"\tNA\tNA\t.";
		last=$4; if ($13!=".") { GID[$26]=$25; TID[$16]=$16; GNA[$27]=$27; no_ol=0;} else no_ol=1;
	} 
	else { GID[$26]=$25; TID[$16]=$16; GNA[$27]=$27; no_ol=0;} 
}' | sed 's/,\t/\t/g' | sed 's/,$//' | 
## we now use awk(above)  instaed groupBy as groupBy uses huge amount of memory for big files!!!
#${bedtools}groupBy -g 4 -c 26,27,16,5,25 -o distinct,distinct,distinct,distinct,mean |
##
## output from "groupBy" is:
##
## SN7001180:281:C99CMACXX:2:2201:15003:5111:SC:CTGTTG:37:UMI:AACCGT:35:25:38	NA	NA	.
## SN7001180:281:C99CMACXX:2:2201:15070:4540:SC:GTGACA:37:UMI:ATCTTA:35:30:38	ENSMUSG00000025902.13	Sox17__chr1	ENSMUST00000192650.5,ENSMUST00000027035.9,ENSMUST00000195555.1
##
## now we can count for each cell and UMI the features, also we make a little QC summary
## we only count a read if it maps to a unique feature/gene, in awk we keep track of the UMIs per gene and output all this for each gene at the end
##
awk -v map_f=$bc_file '
BEGIN{
	while(getline<map_f) {                      ## read in cell barcodes
		CELL[$2]=$1; num_cells+=1;
	}
}
{
	pos=match($1,":SC:");                       ## get startpos of precious barcode info from readname
	split(substr($1,pos+1),BC,":");             ## split on ":" to separate all info and stor in array "BC"
	num=split($2,GENES,",");                    ## get number of features a read overlaps by splitting on ","
	if ( $2!="NA" && num==1 && BC[2] in CELL) { ## if unique feature and a "valid" cell barcode then count it
		ALL[$3][BC[5]][CELL[BC[2]]] += 1;     ## $3 is single GENEID if num==1
		feat_uniq+=1;                         ## only stats
	}
	if (BC[2] in CELL) {                        ## only stats
		if ($2=="NA") cell_nofeat+=1;
			else
				if (num>1) cell_multi+=1;
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
	sum_reads = feat_uniq + cell_nofeat + cell_multi + nocell;
	sum = "FEATURE_UNIQUE\t"feat_uniq"\t"(feat_uniq/sum_reads*100)"\n";
	sum = sum"FEATURE_MULTI\t"cell_multi"\t"(cell_multi/sum_reads*100)"\n";
	sum = sum"CELL_NOFEATURE\t"cell_nofeat"\t"(cell_nofeat/sum_reads*100)"\n";
	sum = sum"NOCELL\t"nocell"\t"(nocell/sum_reads*100)"\n";
	sum = sum"NUM_READS\t"sum_reads"\t100.0";
	print sum > "/dev/stderr";                  ## prints stats to stderr
}' 

