### add barcodes from R1 to R2 #########

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

rule sc_hisat2_genomic:
    input:
        read_barcoded = fastq_dir+"/{sample}"+".fastq.gz"
    output:
        bam = "HISAT2_genomic/{sample}.bam"
    threads: 20
    shell: 
        hisat2_path + " --rna-strandness F -k 5"
        " -x " + hisat2_index + ""
        " -U " + {input.read_barcoded}+ ""
        " --known-splicesite-infile " + known_splicesites + ""
        " --no-unal -p " + {threads} + " --reorder | "
        "grep -P '^@ | NH:i:1\b' | "
        ""+samtools_path + "samtools view -F256 -Sb - | "
        ""+samtools_path + "samtools sort -@5 -m 2G - " {output.bam}