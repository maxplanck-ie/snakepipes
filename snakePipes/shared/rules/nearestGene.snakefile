sample_name = os.path.splitext(os.path.basename(sampleSheet))[0]
if pipeline in ['chip-seq','ATAC-seq']:
    change_direction = ["UP","DOWN","MIXED"]
rule get_nearest_gene:
    input:
        bed=lambda wildcards: "CSAW_{}".format(sample_name)+"/Filtered.results.{change_dir}.bed"  if pipeline in ['chip-seq','ATAC-seq'] else ""
    output:
        annotated_bed="AnnotatedResults_{}".format(sample_name)+"/Filtered.results.{change_dir}_withNearestGene.bed" 
    params:
        genes_bed=genes_bed
    log:
        err= "AnnotatedResults_{}".format(sample_name)+"/logs/bedtools_closest.{change_dir}.err",
    conda: CONDA_RNASEQ_ENV
    shell: "if [ -r {input.bed} ]; then bedtools closest -D b -a <( bedtools sort -i {input.bed} ) -b <( bedtools sort -i {params.genes_bed} ) > {output.annotated_bed};fi 2> {log.err}"
