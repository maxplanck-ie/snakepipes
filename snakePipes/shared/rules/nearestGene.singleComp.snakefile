sample_name = os.path.splitext(os.path.basename(sampleSheet))[0]
if pipeline in ['chipseq','ATACseq']:
    change_direction = ["UP","DOWN","MIXED"]
rule get_nearest_transcript:
    input:
        bed="CSAW_{}_{}".format(peakCaller, sample_name)+"/Filtered.results.{change_dir}.bed"  if pipeline in ['chipseq','ATACseq'] else ""
    output:
        annotated_bed=temp("AnnotatedResults_{}_{}".format(peakCaller, sample_name)+"/Filtered.results.{change_dir}_withNearestTranscript.bed")
    params:
        genes_bed=genes_bed,
        field_offset=lambda wildcards: "18" if pipeline in ['chipseq','ATACseq'] else ""
    conda: CONDA_RNASEQ_ENV
    shell: """
            if [ -r {input.bed} ]; then bedtools closest -D b -a <( bedtools sort -i {input.bed} ) -b <( bedtools sort -i {params.genes_bed} ) | cut -f1-{params.field_offset},$(( {params.field_offset} + 1 ))-$(( {params.field_offset} + 4 )),$(( {params.field_offset} + 6 )),$(( {params.field_offset} + 13 )) > {output.annotated_bed};fi
           """

rule get_nearest_gene:
    input:
        bed="AnnotatedResults_{}_{}".format(peakCaller, sample_name)+"/Filtered.results.{change_dir}_withNearestTranscript.bed",
        t2g="Annotation/genes.filtered.t2g",
        gene_symbol="Annotation/genes.filtered.symbol"
    output:
        annotated_bed="AnnotatedResults_{}_{}".format(peakCaller, sample_name)+"/Filtered.results.{change_dir}_withNearestGene.txt"
    params:
        pipeline=pipeline,
        wdir="AnnotatedResults_{}_{}".format(peakCaller, sample_name)
    conda: CONDA_RNASEQ_ENV
    script: "../rscripts/nearestGene.R"
