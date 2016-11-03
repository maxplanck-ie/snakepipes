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
         temp_dir = os.getenv("TMPDIR", outdir)
     output:
         bam = "HISAT2_genomic/{sample}.bam",
         align_summary = "HISAT2_genomic/{sample}.HISAT2_genomic_summary.txt",
     threads: 20
     shell: 
         "mkdir -p {input.temp_dir} && " + hisat2_path + " --rna-strandness F -k 5"
         " -x " + hisat2_index + ""
         " -U {input.read_barcoded} "
         " --known-splicesite-infile " + known_splicesites + ""
         " --no-unal -p {threads} --reorder 2> {output.align_summary} | "
         "grep -P '^@|NH:i:1\b' | "
         ""+samtools_path + "samtools view -F256 -Sb - | "
         ""+samtools_path + "samtools sort -T {input.temp_dir}/{wildcards.sample} -@5 -m 2G -O bam - > {output.bam}"
         
         
#### count reads/UMIs per gene
#rule make_bed12_from_gtf:
#    input: genes_gtf
#    output: "annotations/genes.bed"
#    run:    if provided_bed != None:
#                shell("cp " + provided_bed + " annotations/genes.bed")
#            else:
                
    
    
# rule sc_get_counts_genomic:
#     input: 
#         bam = "HISAT2_genomic/{sample}.bam",
#         bed12 = "/data/repository/organisms/GRCm38_ensembl/gencode/m9/genes.bed"
#     output: 
#         "Counts/{sample}.counts"
#     shell: 
#         bedtools_path #+ "/intersectBed -a {input.bam} -b <(cat {input.bed12} | grep -v -e "PATCH" -e "CHR" ) -split -bed -wao -s -f 0.9 | awk -v map_f=/data/repository/organisms/GRCm38_ensembl/gencode/m9/transcript2gene.txt 'BEGIN{while (getline < map_f) MAP[$1]=$2}{if ($13!="."){$16=MAP[$16]; OFS="\t";print $0}}' | /package/bedtools2-2.25.0/bin/groupBy -g 4 -c 16,11 -o distinct | less
        