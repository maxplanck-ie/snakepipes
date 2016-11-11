rule filter_exclude_annotation:
    input:
        bed_annot = "Annotation/genes.annotated.bed",
    output:
        bed_filtered = "Annotation/genes.filtered.bed"
    params:
        exclude_pattern =  transcripts_exclude
    shell:
        """ cat {input.bed_annot} | grep -v -P "{params.exclude_pattern}" > {output.bed_filtered}; """
