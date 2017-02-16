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
        read_barcoded = fastq_dir+"/{sample}.fastq.gz",
    output:
        bam = "HISAT2_genomic/{sample}.bam",
        align_summary = "HISAT2_genomic/{sample}.HISAT2_genomic_summary.txt",
    params:
        hisat2_opts = "--pen-cansplice 3 --mp 4,2"
    threads:
        20
    shell: 
        hisat2_path + "hisat2 {params.hisat2_opts} --rna-strandness F -k 5"
        " -x " + hisat2_index + ""
        " -U {input.read_barcoded} "
        " --known-splicesite-infile " + known_splicesites + ""
        " --no-unal -p {threads} --reorder 2> {output.align_summary} | "
        "grep -P '^@|NH:i:1\\b' | "
        ""+samtools_path + "samtools view -F256 -Sb - | "
        ""+samtools_path + "samtools sort -T ${{TMPDIR}}{wildcards.sample} -@5 -m 2G -O bam - > {output.bam}; "
        ""+samtools_path + "samtools index {output.bam} "
        
## STAR genomic mapping 
rule sc_STAR_genomic:
    input:
        read_barcoded = fastq_dir+"/{sample}.fastq.gz",
        gtf = "Annotation/genes.filtered.gtf"
    output:
        bam = "STAR_genomic/{sample}.bam"
    params:
        opts = "--sjdbOverhang 100 --twopassMode Basic"
    threads:
        20
    shell:
        star_path + "STAR --genomeDir "+star_index + ""
        " --runThreadN {threads} --readFilesIn {input.read_barcoded} "
        " --readFilesCommand zcat --outFileNamePrefix ${{TMPDIR}}/{wildcards.sample}. " 
        " --sjdbGTFfile {input.gtf} {params.opts} --outStd SAM | "
        #" grep -P '^@|NH:i:1\\b' | "
        "" + samtools_path + "samtools sort -T ${{TMPDIR}}tmp_{wildcards.sample} -@5 -m 2G -O bam - > {output.bam}; "
        "" + samtools_path + "samtools index {output.bam}; "
        " cp ${{TMPDIR}}/{wildcards.sample}.Log.final.out STAR_genomic/;"


rule sc_bam_featureCounts_genomic:
    input:
        bam = "STAR_genomic/{sample}.bam",
        gtf = "Annotation/genes.filtered.gtf"
    output:
        counts = "Counts/{sample}.cout.csv",
        counts_summary = "Counts/{sample}.cout_summary.txt",
        cell_summary = "Counts/{sample}.cout_summary.cells.txt"
    params:
        count_script = workflow.basedir+"/scRNAseq_bam_featureCounts.sh",
        bc_file = barcode_file,
        fc_path = feature_counts_path
    threads: 
        5
    shell:
        """
        {params.count_script} {input.bam} {input.gtf} {params.bc_file} {wildcards.sample} {params.fc_path} ${{TMPDIR}} {threads} 1>{output.counts} 2>{output.counts_summary};
        cat {output.counts_summary} | sed -n -e '/#idx/,$p' > {output.cell_summary} 
        """

# rule sc_get_counts_genomic:
#     input:
#         bam = "STAR_genomic/{sample}.bam",
#         bed = "Annotation/genes.filtered.bed"
#     output: 
#         counts = "Counts/{sample}.cout.csv",
#         counts_summary = "Counts/{sample}.cout_summary.txt"
#     params:
#         count_script = workflow.basedir+"/scRNAseq_bam_genomic_feature_count.sh",
#         bc_file = barcode_file,
#         bedtools = bedtools_path,
#         samtools = samtools_path
#     threads:
#         5
#     shell: 
#         """
#             {params.count_script} {input.bam} {input.bed} {params.bc_file} {params.bedtools} {params.samtools} 1>{output.counts} 2>{output.counts_summary}
#         """


rule scale_counts:
    input:
        counts = "Counts/{sample}.cout.csv"
    output:
        coutt = "Counts/{sample}.coutt.csv",
        coutb = "Counts/{sample}.coutb.csv",
        coutc = "Counts/{sample}.coutc.csv",
        log = "Counts/{sample}.extract_log.txt",
    params:
        count_script = workflow.basedir+"/extract_counts_rb.pl",
        UMI_length = UMI_length
    shell:
        """{params.count_script} -bl={params.UMI_length} -in={input.counts} """
        """ -outc={output.coutc} -outb={output.coutb} -outt={output.coutt} &>{output.log} """


rule combine_sample_counts:
    input:
        expand("Counts/{sample}.coutt.csv",sample = samples)
    output:
        merged_matrix = "Results/all_samples.gencode_genomic.coutt_merged.csv"
    params:
        merge_script = workflow.basedir+"/scRNAseq_merge_coutt_files2.R",
        split = 1
    shell:
        R_path+"""Rscript {params.merge_script} Counts/ {output.merged_matrix} {params.split} """


rule sc_QC_metrics:
    input:
        expand("Counts/{sample}.coutb.csv",sample = samples),
        expand("Counts/{sample}.coutc.csv",sample = samples)
    output:
        summary = "QC_report/QC_report.all_samples.tsv",
        summary_nice = "QC_report/QC_report.all_samples.txt",
        sc_dat = "QC_report/all_samples.cellsum_coutc_coutb.tsv"
    params: 
        in_dir = outdir+"/Counts/",
        out_dir = outdir+"/QC_report/",
        plot_script = workflow.basedir+"/scRNAseq_QC_metrics.R"
    shell:
        ""+workflow.basedir+"/scRNAseq_QC_metrics.sh {params.in_dir} {params.out_dir} 1>{output.summary};"
        "cat {output.summary} | column -t > {output.summary_nice};"
        ""+R_path+"Rscript {params.plot_script} {output.sc_dat};"


rule bamCoverage_RPKM:
    input:
        bam = "STAR_genomic/{sample}.bam"
    output:
        "Tracks/{sample}.Coverage.bw"
    params:
        bw_binsize = config["bw_binsize"]
    log:
        "Tracks/logs/bamCoverage_Coverage.{sample}.log"
    benchmark:
        "Tracks/.benchmark/bamCoverage_Coverage.{sample}.benchmark"
    threads: 8
    shell:
        deepTools_path+"bamCoverage "
        "-b {input.bam} "
        "-o {output} "
        "--binSize {params.bw_binsize} "
        "-p {threads} "
        "&> {log}"


rule plotEnrichment:
    input:
        bam = expand("STAR_genomic/{sample}.bam", sample=samples),
        bed = "Annotation/genes.filtered.bed"
    output:
        png = "deepTools_qc/plotEnrichment/plotEnrichment.png",
        tsv = "deepTools_qc/plotEnrichment/plotEnrichment.tsv",
    params:
        labels = " ".join(samples),
    log:
        "deepTools_qc/logs/plotEnrichment.log"
    benchmark:
        "deepTools_qc/.benchmark/plotEnrichment.benchmark"
    threads: 8
    shell:
        deepTools_path+"plotEnrichment "
        "-p {threads} "
        "-b {input.bam} "
        "--BED {input.bed} "
        "--plotFile {output.png} "
        "--labels {params.labels} "
        "--plotTitle 'Fraction of reads in regions' "
        "--outRawCounts {output.tsv} "
        "--variableScales "
        "&> {log} "

## zcat 14wks_Eed_WT_1.umi.fastq.gz | /package/hisat2-2.0.4/hisat2 --rna-strandness F -k 5 -x /data/repository/organisms/GRCm38_ensembl/HISAT2Index/genome -U - --no-unal -p 16 --reorder | grep -P '^@|NH:i:1\b' | samtools view -F256 -Sb - | /package/bedtools2-2.25.0/bin/intersectBed -a - -b <(cat /data/repository/organisms/GRCm38_ensembl/gencode/m9/genes.bed| grep -v -e "PATCH" -e "CHR" ) -split -bed  -wo -s | awk -v map_f=gencode.M9.full.table 'BEGIN{while (getline < map_f) {MAP[$2]=$1;MAP2[$2]=$4}}{if ($13!="."){OFS="\t";print $0,MAP[$16],MAP2[$16]"__chr"$1}}' | /package/bedtools2-2.25.0/bin/groupBy -g 4 -c 26,27,16,5,25 -o distinct,distinct,distinct,collapse,mean | awk -v map_f=/data/pospisilik/group/heyne/scRNAseq/sagar/celseq_barcodes.192.txt 'BEGIN{while (getline < map_f) {CELL[$2]=$1;COUNTS[$1]=0}}{pos=match($1,":SC:");split(substr($1,pos+1),BC,":"); num=split($2,GENES,",");if ( (num==1 && BC[2] in CELL) ) {if (!($3 in ALL)){for (i=1; i<=192;i++) ALL[$3][i]=0;} ALL[$3][CELL[BC[2]]] += 1}}END{for (i in ALL){printf i" "; for (j=1;j<=192;j++){ printf ALL[i][j]" ";} printf "\n"}}' | less
##     cat test.bam | /package/bedtools2-2.25.0/bin/intersectBed -a - -b <(cat /data/repository/organisms/GRCm38_ensembl/gencode/m9/genes.bed | grep -v -e "PATCH" -e "CHR" ) -split -bed  -wo -s | awk -v map_f=gencode.M9.full.table 'BEGIN{while (getline < map_f) {MAP[$2]=$1;MAP2[$2]=$4}}{if ($13!="."){OFS="\t";print $0,MAP[$16],MAP2[$16]"__chr"$1}}' | /package/bedtools2-2.25.0/bin/groupBy -g 4 -c 26,27,16,5,25 -o distinct,distinct,distinct,collapse,mean | awk -v map_f=/data/pospisilik/group/heyne/scRNAseq/sagar/celseq_barcodes.192.txt 'BEGIN{while (getline < map_f) {CELL[$2]=$1;COUNTS[$1]=0}}{pos=match($1,":SC:");split(substr($1,pos+1),BC,":"); num=split($2,GENES,",");if ( (num==1 && BC[2] in CELL) ) {if (!($3 in ALL) || !(BC[5] in ALL[$3])){for (i=1; i<=192;i++) ALL[$3][BC[5]][i]=0;} ALL[$3][BC[5]][CELL[BC[2]]] += 1}}END{for (i in ALL){for (k in ALL[i]){printf i" "k" "; for (j=1;j<=192;j++){ printf ALL[i][k][j]" ";} printf "\n"}}}' > test.cout.csv
#bedtools_path #+ "/intersectBed -a {input.bam} -b <(cat {input.bed12} | grep -v -e "PATCH" -e "CHR" ) -split -bed -wao -s -f 0.9 | awk -v map_f=/data/repository/organisms/GRCm38_ensembl/gencode/m9/transcript2gene.txt 'BEGIN{while (getline < map_f) MAP[$1]=$2}{if ($13!="."){$16=MAP[$16]; OFS="\t";print $0}}' | /package/bedtools2-2.25.0/bin/groupBy -g 4 -c 16,11 -o distinct | less
