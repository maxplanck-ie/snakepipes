## make statdard annotation

if genome.startswith(("hg","mm")):
    rule create_annotation_bed:
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
            """ OFS="\\t"; print tid,tna,tt,gid,gna,gt,"gencode",basic,tsl,lvl}}' | """
            """ tr " " "\\t" | sort -k2) | """
            """ awk '{{$13=$13"\\t"$1; $4=$4"\\t"$1; OFS="\\t";print $0}}' | """
            """ cut --complement -f 1,14,16,18,20,22,24 > {output.bed_annot} """
elif genome.startswith(("dm")):
    rule create_annotation_bed:
        input:
            gtf = genes_gtf
        output:
            bed_annot = "Annotation/genes.annotated.bed"
        shell:
            "join -t $'\t' -o auto --check-order -1 4 -2 2 "
            "<("+UCSC_tools_path+"gtfToGenePred {input.gtf} /dev/stdout | "+UCSC_tools_path+"""genePredToBed /dev/stdin /dev/stdout | tr " " "\\t" | sort -k4) """
            """ <(cat {input.gtf} | awk '$3=="transcript"{{print $0}}' | tr -d "\\";" | """
            """ awk '{{"""
            """ pos=match($0,"gene_biotype.[^[:space:]]+"); gt=substr($0,RSTART,RLENGTH); """
            """ pos=match($0,"transcript_biotype.[^[:space:]]+");tt=substr($0,RSTART,RLENGTH); """
            """ pos=match($0,"gene_id.[^[:space:]]*"); gid=substr($0,RSTART,RLENGTH); """
            """ pos=match($0,"transcript_id.[^[:space:]]*"); tid=substr($0,RSTART,RLENGTH); """
            """ pos=match($0,"transcript_name.[^[:space:]]*"); tna=substr($0,RSTART,RLENGTH); """
            """ pos=match($0,"gene_name.[^[:space:]]*"); gna=substr($0,RSTART,RLENGTH); """
            """ OFS="\\t"; print tid,tna,tt,gid,gna,gt}}' | """
            """ tr " " "\\t" | sort -k2) | """
            """ awk '{{$13=$13"\\t"$1; $4=$4"\\t"$1; OFS="\\t";print $0}}' | """
            """ cut --complement -f 1,14,16,18,20,22,24 > {output.bed_annot} """


rule filter_annotation_bed:
    input:
        bed_annot = "Annotation/genes.annotated.bed",
    output:
        bed_filtered = "Annotation/genes.filtered.bed"
    params:
        pattern =  filter_annotation
    shell:
        "cat {input.bed_annot} | grep {params.pattern} > {output.bed_filtered} "
