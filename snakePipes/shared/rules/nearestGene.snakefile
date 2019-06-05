sample_name = os.path.splitext(os.path.basename(sampleSheet))[0]
change_direction = ["UP","DOWN","MIXED"]
rule get_nearest_gene:
    input:
        bed=lambda wildcards: "CSAW_{}/Filtered.results.{change_dir}.bed".format(sample_name) if pipeline in ['chip-seq','ATAC-seq']
    output:
        annotated_bed=lambda wildcards: "CSAW_{}/Filtered.results.{change_dir}_withNearestGene.bed".format(sample_name) if pipeline in ['chip-seq','ATAC-seq']
    params:
        out_log=lambda wildcards: "CSAW_{}/logs/bedtools_closest.{change_dir}.out".format(sample_name) if pipeline in ['chip-seq','ATAC-seq'],
        err_log=lambda wildcards: "CSAW_{}/logs/bedtools_closest.{change_dir}.err".format(sample_name) if pipeline in ['chip-seq','ATAC-seq'],
        genes_bed=genes_bed
    conda: CONDA_RNASEQ_ENV
    shell: "if [ -r {input.bed} ]; then bedtools closest -D b -a <( sort -k1,1 -k2,2n {input.bed} ) -b <( sort -k1,1 -k2,2n {params.genes_bed} ) > {output.annotated_bed};fi 2> {params.err_log}"
