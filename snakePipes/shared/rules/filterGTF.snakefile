## make standard annotation

## rule creates a bed12+ file based on the UCSC tools
## we also add 12 annotation columns that are helpful to get gene names, biotypes etc. from gene ids or for filtering
##     cols 13-15:    transcript_id, transcript_name, transcript_biotype
##     cols 16-18:    gene_id, gene_name, gene_biotype
##
## and some gencode specific tags if present:
##     cols 19+20:    gencode [basic, full, NA], basic=representative transcripts for each gene
##     cols 21+22:    transcript_support_level [NA,1,2,3,4,5] 
##        1 (all splice junctions of the transcript are supported by at least one non-suspect mRNA),
##        2 (the best supporting mRNA is flagged as suspect or the support is from multiple ESTs),
##        3 (the only support is from a single EST),
##        4 (the best supporting EST is flagged as suspect),
##        5 (no single transcript supports the model structure),
##        NA (the transcript was not analyzed)  
##     cols 23+24:    level [1,2,3,NA]  1=verified loci, 2=manually annotated loci, 3=automatically annotated loci
##
## according to the gtf format description, gene_id/transcript_id tags are required (by UCSCStools)
## missing tags we set to NA or to the gene_id/transcript_id 
rule create_annotation_bed:
    input:
        gtf = genes_gtf
    output:
        bed_annot = "Annotation/genes.annotated.bed"
    conda: CONDA_RNASEQ_ENV
    shell:
        "join -t $'\t' -o auto --check-order -1 4 -2 1 "
        "<(gtfToGenePred -ignoreGroupsWithoutExons {input.gtf} /dev/stdout | genePredToBed /dev/stdin /dev/stdout | tr ' ' $'\t' | sort -k4,4) "
        "<(cat {input.gtf} | "
        """ awk '$1!~"^#" && $3~"transcript|exon|start_codon|stop_codon|CDS"{{print $0}}' | tr -d "\\";" | """
        """ awk '{{pos=match($0,"gene_id.([^[:space:]]+)",a); if (pos!=0) gid=a[1];else exit 1; """
        """ pos=match($0,"transcript_id.([^[:space:]]+)",a); if (pos!=0) tid=a[1]; else tid=gid; """
        """ pos=match($0,"transcript_name.([^[:space:]]+)",a); if (pos!=0) tna=a[1]; else tna=tid; """
        """ pos=match($0,"transcript_[bio]*type.([^[:space:]]+)",a); if (pos!=0) tt=a[1]; else tt="NA"; """
        """ pos=match($0,"gene_name.([^[:space:]]+)",a);  if (pos!=0) gna=a[1]; else gna=gid; """
        """ pos=match($0,"gene_[bio]*type.([^[:space:]]+)",a); if (pos!=0) gt=a[1]; else gt="NA"; """
        """ pos=match($0,"transcript_support_level.([^[:space:]]+)",a); if (pos!=0) tsl=a[1]; else tsl="NA"; """
        """ pos=match($0,"[[:space:]]level.([^[:space:]]+)",a); if (pos!=0) lvl=a[1] ; else lvl="NA"; """
        """ pos=match($0,"tag.basic"); if (lvl!~"NA"){{if (pos==0) basic="full"; else basic="basic"}} else basic="NA"; """
        """ OFS="\\t"; print tid,tna,tt,gid,gna,gt,"gencode",basic,"transcript_support_level",tsl,"level",lvl}}' | """
        """ sort | uniq | sort -k1,1) | """
        """ awk '{{OFS="\\t";print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$1,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}}' > {output.bed_annot} """ 


## filter transcripts with grep
rule filter_annotation_bed:
    input:
        bed_annot = "Annotation/genes.annotated.bed"
    output:
        bed_filtered = "Annotation/genes.filtered.bed"
    params:
        pattern = str(filterGTF or '\'\''),
        cmd = "Annotation/filter_command.txt"
    shell:
        """
        cat {input.bed_annot} | grep {params.pattern} > {output.bed_filtered};
        echo 'cat {input.bed_annot} | grep \'{params.pattern}\' > {output.bed_filtered}' > {params.cmd}
        """

## t2g files is used in alignment-free mode for sleuth/DESeq2
## transcript_id -> gene_id mapping table
rule annotation_bed2t2g:
    input:
        bed_annot = 'Annotation/genes.filtered.bed' 
    output:
        'Annotation/genes.filtered.t2g'
    shell:
        "cat {input.bed_annot} | awk '{{OFS=\"\t\"; print $13,$16,$14}}' > {output}"


## used in DESeq2, gene_id -> gene_name mapping table
rule annotation_bed2symbol:
    input:
        bed_annot = 'Annotation/genes.filtered.bed' 
    output:
        'Annotation/genes.filtered.symbol'
    shell:
        "cat {input.bed_annot} | awk '{{OFS=\"\t\"; print $16, $17}}' | sort | uniq > {output}"


## creates fasta sequence file for the filtered transcripts
## used to create the salmon index in alignment-free mode 
rule annotation_bed2fasta:
    input:
        bed = "Annotation/genes.filtered.bed",
        genome_fasta = genome_fasta
    output:
        "Annotation/genes.filtered.fa"
    benchmark:
        "Annotation/.benchmark/annotation_bed2fasta.benchmark"
    threads: 1
    conda: CONDA_RNASEQ_ENV
    shell:
        "bedtools getfasta -name -s -split -fi {input.genome_fasta} -bed <(cat {input.bed} | cut -f1-12) | sed 's/(.*)//g' > {output}"


## the gtf created here is used for featureCounts, contains only the filtered entries from original gtf  
rule annotation_bed2gtf:
    input:
        bed = "Annotation/genes.filtered.bed"
    output:
        gtf = "Annotation/genes.filtered.gtf"
    conda: CONDA_RNASEQ_ENV
    shell:
    	"""
        bedToGenePred {input.bed} stdout | awk -v map_f={input.bed} '
        BEGIN{{while (getline < map_f) MAP[$13]=$16}} {{OFS="\\t";print $0,"0",MAP[$1]}}' |
        genePredToGtf -utr file stdin stdout |
        awk -v map_f={input.bed} '
        BEGIN{{while (getline < map_f) MAP[$16]=$17}}
        {{match($0,"gene_name[[:space:]\\";]+([^[:space:]\\";]+)",a);
        sub("gene_name[[:space:];\\"]+"a[1],"gene_name \\""MAP[a[1]],$0); print $0}}' > {output.gtf} 
        """ 
