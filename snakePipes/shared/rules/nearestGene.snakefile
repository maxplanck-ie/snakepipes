sample_name = os.path.splitext(os.path.basename(sampleSheet))[0]
if pipeline in ['chip-seq','ATAC-seq']:
    change_direction = ["UP","DOWN","MIXED"]
rule get_nearest_transcript:
    input:
        bed="CSAW_{}_{}".format(peakCaller, sample_name)+"/Filtered.results.{change_dir}.bed"  if pipeline in ['chip-seq','ATAC-seq'] else ""
    output:
        annotated_bed=temp("AnnotatedResults_{}_{}".format(peakCaller, sample_name)+"/Filtered.results.{change_dir}_withNearestTranscript.bed")
    params:
        genes_bed=genes_bed,
        field_offset=lambda wildcards: "18" if pipeline in ['chip-seq','ATAC-seq'] else ""
    log:
        err= "AnnotatedResults_{}_{}".format(peakCaller, sample_name)+"/logs/bedtools_closest.{change_dir}.err",
    conda: CONDA_RNASEQ_ENV
    shell: """
            if [ -r {input.bed} ]; then bedtools closest -D b -a <( sort -k 1,1 -k2,2n {input.bed} | tr \" \" \"\t\"  ) -b <( sort -k 1,1 -k2,2n {params.genes_bed} | tr \" \" \"\t\" ) | cut -f1-{params.field_offset},$(( {params.field_offset} + 1 ))-$(( {params.field_offset} + 4 )),$(( {params.field_offset} + 6 )),$(( {params.field_offset} + 13 )) > {output.annotated_bed};fi 2> {log.err}
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
    log:
        err="AnnotatedResults_{}_{}".format(peakCaller, sample_name)+"/logs/nearestGene.{change_dir}.err",
        out="AnnotatedResults_{}_{}".format(peakCaller, sample_name)+"/logs/nearestGene.{change_dir}.out"
    conda: CONDA_RNASEQ_ENV
    script: "../rscripts/nearestGene.R"
