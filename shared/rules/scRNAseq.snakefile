### add barcodes from R1 to R2 #########
#from xdiagnose.utils.processes import shell

rule fastq_barcode:
        input:
            R1 = "FASTQ/{sample}"+reads[0]+".fastq.gz",
            R2 = "FASTQ/{sample}"+reads[1]+".fastq.gz"
        output:
            R2_barcoded = "FASTQ_barcoded/{sample}.fastq.gz"
        params:
            UMI_length = UMI_length,
            UMI_offset = UMI_offset,
            CELLI_length = CELLI_length,
            CELLI_offset = CELLI_offset
        benchmark:
            "FASTQ_barcoded/.benchmark/fastq_barcoded.{sample}.benchmark"
        threads: 2
        shell:"""
            paste <(paste - - - - < <(zcat {input.R1}))   <(paste - - - - < <(zcat {input.R2})) | \
            tr '\t' '\n' | \
            awk -v CBAR_LEN={params.CELLI_length} -v UMI_LEN={params.UMI_length} \
                -v CBAR_OFFSET={params.CELLI_offset} -v UMI_OFFSET={params.UMI_offset} ' \
           	BEGIN{{
				for(n=0;n<256;n++) 
					phred33[sprintf("%c",n)]=n-33
 			}}
			{{
			if (NR%8==2) 
		 		{{CB=substr($0,CBAR_OFFSET,CBAR_LEN);
	 			UMI=substr($0,UMI_OFFSET,UMI_LEN);
	 			if (CBAR_OFFSET+CBAR_LEN>UMI_OFFSET+UMI_LEN)
	 			TAIL=substr($0,CBAR_OFFSET+CBAR_LEN); else
	 			TAIL=substr($0,UMI_OFFSET+UMI_LEN);
	 			NUMT=gsub(/T/,"#",TAIL);
			}} 
			if (NR%8==4) 
				{{split(substr($0,UMI_LEN+1,CBAR_LEN),CBQA,"");
	 			split(substr($0,1,UMI_LEN),UMIQA,""); 
	 			QUAL_UMI=0; QUAL_CB=0;
 	 			for (i=1;i<=length(UMIQA);i++)
					{{QUAL_UMI+=phred33[UMIQA[i]]}}; 
	 			for (i=1;i<=length(CBQA);i++)
					{{QUAL_CB+=phred33[CBQA[i]]}};
				}}; 
			OFS=" ";
			if (NR%8==5) 
				{{$1=$1":SC:"CB":"sprintf("%.0f",QUAL_CB/CBAR_LEN)":UMI:"UMI":"sprintf("%.0f",QUAL_UMI/UMI_LEN)":"NUMT":"length(TAIL);print $0;
				}}; 
			if (NR%8==0 || NR%8>5) print $0}}' | pigz -c -p 2 > {output.R2_barcoded}
            """

### HISAT2 genomic mapping
rule sc_hisat2_genomic:
     input:
         read_barcoded = fastq_dir+"/{sample}"+".fastq.gz",
     output:
         bam = "HISAT2_genomic/{sample}.bam",
         align_summary = "HISAT2_genomic/{sample}.HISAT2_genomic_summary.txt",
     threads: 20
     shell: 
         hisat2_path + "hisat2 --rna-strandness F -k 5"
         " -x " + hisat2_index + ""
         " -U {input.read_barcoded} "
         " --known-splicesite-infile " + known_splicesites + ""
         " --no-unal -p {threads} --reorder 2> {output.align_summary} | "
         "grep -P '^@|NH:i:1\\b' | "
         ""+samtools_path + "samtools view -F256 -Sb - | "
         ""+samtools_path + "samtools sort -T ${{TMPDIR}}{wildcards.sample} -@5 -m 2G -O bam - > {output.bam}"
         
         
#### count reads/UMIs per gene
#rule make_bed12_from_gtf:
#    input: genes_gtf
#    output: "annotations/genes.bed"
#    run:    if provided_bed != None:
#                shell("cp " + provided_bed + " annotations/genes.bed")
#            else:
                
## make statdard annotation
rule create_annotation:
    input: 
        gtf = genes_gtf
    output: 
        bed_annot = "Annotation/genes.annotated.bed"
    shell:
        "join -t $'\t' -o auto --check-order -1 4 -2 2 "
        "<("+UCSC_tools_path+"gtfToGenePred {input.gtf} /dev/stdout | "+UCSC_tools_path+"""genePredToBed /dev/stdin /dev/stdout | tr " " "\\t" | sort -k4) """
         """ <(cat {input.gtf} | awk '$3=="transcript"{{print $0}}' | tr -d "\\";" | """
         """ awk '{{pos=match($0,"tag.basic"); if (pos==0) basic="full"; else basic="basic"; """
         """ pos=match($0,"gene_type.[^[:space:]]+"); gt=substr($0,RSTART,RLENGTH); """
         """ pos=match($0,"transcript_type.[^[:space:]]+");tt=substr($0,RSTART,RLENGTH); """
         """ pos=match($0,"transcript_support_level.[^[:space:]]+"); if (pos!=0) tsl=substr($0,RSTART,RLENGTH);else tsl="transcript_support_level NA"; """
         """ pos=match($0,"[[:space:]]level.[^[:space:]]*"); lvl=substr($0,RSTART,RLENGTH); """
         """ pos=match($0,"gene_id.[^[:space:]]*"); gid=substr($0,RSTART,RLENGTH); """
         """ pos=match($0,"transcript_id.[^[:space:]]*"); tid=substr($0,RSTART,RLENGTH); """
         """ pos=match($0,"transcript_name.[^[:space:]]*"); tna=substr($0,RSTART,RLENGTH); """
         """ pos=match($0,"gene_name.[^[:space:]]*"); gna=substr($0,RSTART,RLENGTH); """
         """ OFS="\\t"; print tid,tna,gid,gna,"gencode",basic,tt,gt,tsl,lvl}}' | """
         """ tr " " "\\t" | sort -k2) | """
         """ awk '{{$13=$13"\\t"$1; $4=$4"\\t"$1; OFS="\\t";print $0}}' | """
         """ cut --complement -f 1,14,16,18,20 > {output.bed_annot} """
         
rule filter_annotation:
    input:
        bed_annot = "Annotation/genes.annotated.bed",    
    output:
        bed_filtered = "Annotation/genes.filtered.bed"
    params: 
        exclude_pattern =  transcripts_exclude  
    shell:
        """ cat {input.bed_annot} | grep -v -P "{params.exclude_pattern}" > {output.bed_filtered}; """ 
        
rule sc_get_counts_genomic:
    input:
        bam = "HISAT2_genomic/{sample}.bam",
        bed = "Annotation/genes.filtered.bed"
    output: 
        counts = "Counts/{sample}.rbar.counts",
        counts_summary = "Counts/{sample}.rbar.counts_summary.txt"
    params:
        map_tab = "Annotation/genes.filtered.bed",
        bc_file = barcode_file
    shell: 
        bedtools_path+
        """intersectBed -a {input.bam} -b <(cat {input.bed} | cut -f1-12) """
        """ -split -bed -wao -s -nonamecheck | """ 
        """ awk -v map_f={input.bed} """
        """  'BEGIN{{while (getline < map_f) {{MAP[$13]=$15;MAP2[$13]=$16}} }}"""
        """  {{OFS="\\t";if ($13!=".") print $0,MAP[$16],MAP2[$16]"__chr"$1; else print $0,"NA","NA"; }} ' | """
#        """ END{{sum="ASSIGNED:\\t"ass"\\nNOT ASSIGNED:\\t"notass"\\t"(notass/(ass+notass))"\\nALL:\\t"ass+notass; """
#        """  print sum >"/dev/stderr"}} ' 2>{output.counts_summary} | """
        +bedtools_path+
        """groupBy -g 4 -c 26,27,16,5,25 -o distinct,distinct,distinct,collapse,mean | """
        """ awk -v map_f={params.bc_file} """
        """ 'BEGIN{{while(getline<map_f) {{CELL[$2]=$1; COUNTS[$1]=0 }};}} """
        """ {{pos=match($1,":SC:"); split(substr($1,pos+1),BC,":"); num=split($2,GENES,","); """
        """  if ( $2!="NA" && num==1 && BC[2] in CELL) {{ """
        """   if (!($3 in ALL) || !(BC[5] in ALL[$3])) {{ """
        """    for (i=1; i<=192;i++) ALL[$3][BC[5]][i]=0;}} """
        """   ALL[$3][BC[5]][CELL[BC[2]]] += 1; feat_uniq+=1 """
        """  }} """
        """ if (num>1) feat_multi+=1; """ 
        """ if ($2!="NA"){{ """
        """  if (BC[2] in CELL) feat_cell+=1; else feat_nocell+=1;"""
        """ }} else {{"""
        """  if (BC[2] in CELL) nofeat_cell+=1; else nofeat_nocell+=1;}} """
        """ }} """
        """ END{{ """
        """  for (i in ALL) {{ """
        """   for (k in ALL[i]) {{ """
        """     printf i" "k" "; """
        """     for (j=1;j<=192;j++) {{printf ALL[i][k][j]" ";}} """
        """     printf "\\n" """
        """   }} """
        """  }} """
        """ sum_reads = feat_cell + feat_nocell + nofeat_cell + nofeat_nocell; """
        """ sum  = "FEATURE_UNIQUE\\t"feat_uniq"\\t"(feat_uniq/sum_reads)"\\n"; """
        """ sum += "FEATURE_MULTI\\t"feat_multi"\\t"(feat_multi/sum_reads)"\\n"; """
        """ sum += "FEATURE_CELL\\t"feat_cell"\\t"(feat_cell/sum_reads)"\\n"; """
        """ sum += "FEATURE_NOCELL\\t"feat_nocell"\\t"(feat_nocell/sum_reads)"\\n"; """
        """ sum += "NOFEATURE_CELL\\t"nofeat_cell"\\t"(nofeat_cell/sum_reads)"\\n"; """
        """ sum += "NOFEATURE_NOCELL\\t"nofeat_nocell"\\t"(nofeat_nocell/sum_reads)"\\n"; """
        """ sum += "NUM_READS\\t"sum_reads; """
        """ print sum > "/dev/stderr" """
        """ }}' 2>{output.counts_summary} 1> {output.counts} """

## zcat 14wks_Eed_WT_1.umi.fastq.gz | /package/hisat2-2.0.4/hisat2 --rna-strandness F -k 5 -x /data/repository/organisms/GRCm38_ensembl/HISAT2Index/genome -U - --no-unal -p 16 --reorder | grep -P '^@|NH:i:1\b' | samtools view -F256 -Sb - | /package/bedtools2-2.25.0/bin/intersectBed -a - -b <(cat /data/repository/organisms/GRCm38_ensembl/gencode/m9/genes.bed| grep -v -e "PATCH" -e "CHR" ) -split -bed  -wo -s | awk -v map_f=gencode.M9.full.table 'BEGIN{while (getline < map_f) {MAP[$2]=$1;MAP2[$2]=$4}}{if ($13!="."){OFS="\t";print $0,MAP[$16],MAP2[$16]"__chr"$1}}' | /package/bedtools2-2.25.0/bin/groupBy -g 4 -c 26,27,16,5,25 -o distinct,distinct,distinct,collapse,mean | awk -v map_f=/data/pospisilik/group/heyne/scRNAseq/sagar/celseq_barcodes.192.txt 'BEGIN{while (getline < map_f) {CELL[$2]=$1;COUNTS[$1]=0}}{pos=match($1,":SC:");split(substr($1,pos+1),BC,":"); num=split($2,GENES,",");if ( (num==1 && BC[2] in CELL) ) {if (!($3 in ALL)){for (i=1; i<=192;i++) ALL[$3][i]=0;} ALL[$3][CELL[BC[2]]] += 1}}END{for (i in ALL){printf i" "; for (j=1;j<=192;j++){ printf ALL[i][j]" ";} printf "\n"}}' | less
##     cat test.bam | /package/bedtools2-2.25.0/bin/intersectBed -a - -b <(cat /data/repository/organisms/GRCm38_ensembl/gencode/m9/genes.bed | grep -v -e "PATCH" -e "CHR" ) -split -bed  -wo -s | awk -v map_f=gencode.M9.full.table 'BEGIN{while (getline < map_f) {MAP[$2]=$1;MAP2[$2]=$4}}{if ($13!="."){OFS="\t";print $0,MAP[$16],MAP2[$16]"__chr"$1}}' | /package/bedtools2-2.25.0/bin/groupBy -g 4 -c 26,27,16,5,25 -o distinct,distinct,distinct,collapse,mean | awk -v map_f=/data/pospisilik/group/heyne/scRNAseq/sagar/celseq_barcodes.192.txt 'BEGIN{while (getline < map_f) {CELL[$2]=$1;COUNTS[$1]=0}}{pos=match($1,":SC:");split(substr($1,pos+1),BC,":"); num=split($2,GENES,",");if ( (num==1 && BC[2] in CELL) ) {if (!($3 in ALL) || !(BC[5] in ALL[$3])){for (i=1; i<=192;i++) ALL[$3][BC[5]][i]=0;} ALL[$3][BC[5]][CELL[BC[2]]] += 1}}END{for (i in ALL){for (k in ALL[i]){printf i" "k" "; for (j=1;j<=192;j++){ printf ALL[i][k][j]" ";} printf "\n"}}}' > test.cout.csv
#bedtools_path #+ "/intersectBed -a {input.bam} -b <(cat {input.bed12} | grep -v -e "PATCH" -e "CHR" ) -split -bed -wao -s -f 0.9 | awk -v map_f=/data/repository/organisms/GRCm38_ensembl/gencode/m9/transcript2gene.txt 'BEGIN{while (getline < map_f) MAP[$1]=$2}{if ($13!="."){$16=MAP[$16]; OFS="\t";print $0}}' | /package/bedtools2-2.25.0/bin/groupBy -g 4 -c 16,11 -o distinct | less
