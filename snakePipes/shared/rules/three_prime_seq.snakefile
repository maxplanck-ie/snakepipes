# adapted from Andrew Rezansoff's 3' seq pipeline tools
# /data/hilgers/group/rezansoff/3seq_pipeline_tools/

from pathlib import Path
import pandas as pd

# prevent wildcards from picking up directories
wildcard_constraints: sample="[^(\/)]+"

# read in sampleSheet metadata to merge replicates for preprocess_cluster_pas 
sample_metadata = pd.read_table(sampleSheet, index_col=None)

# given a condition, return a list of all samples associated with it
# check that there are not overlapping samples
def samples_from_condition(condition):
    sub_df = sample_metadata[sample_metadata['condition'] == condition]
    assert(set(sub_df['name']) == set(sub_df['name'].unique()))
    return list(sub_df['name'])

# the opposite of the above
def condition_from_sample(sample):
    sub_df = sample_metadata[sample_metadata['name'] == sample]
    return sub_df['condition'].iloc[0]

def get_outdir(folder_name,sampleSheet):
    sample_name = os.path.splitext(os.path.basename(str(sampleSheet)))[0]
    return("{}_{}".format(folder_name, sample_name))


# TODO: 
# 3. check why cmatrix_filtered is commented out - OK
# 5. re-implement clusterPAS?  
# 6. implement filtering of column 4 of the geneAssociation output - things with multiple hits (i.e. commas) should be filtered out [x]
# possibly assign cluster to "next" gene annotation (i.e. the closest?) 
# this would be in signal2gene.py
# 7. assign unique IDs with _0, _1 to 4th column in clusterPAS output in 5'->3' order to create unique 
# example is /data/hilgers/group2/rezansoff/sakshiProj/2507_3seq/polyA_annotation_polysome2K/AllSamples_clustered_strict1_fromSortedNumeric.txt -> 
# /data/hilgers/group2/rezansoff/sakshiProj/2507_3seq/polyA_annotation_polysome2K/AllSamples_clustered_strict1_fromSortedNumeric_uniqIDs.txt
# this will be run by countReadEnds
# then possibly DESeq on countReadEnds output
# understand cmatrix steps? 

tools_dir = Path(maindir) / "shared" / "tools"

# trimming done using STAR and fastp

# rule clusterPAS:
#     input: 
#         "three_prime_seq/SampleAll.txt"
#     conda:
#         CONDA_SHARED_ENV
    

rule polyAT: 
    input: 
        two_bit=genome_2bit,
        bed=genes_bed
    output: 
        "three_prime_seq/poly{base}.bed"
    params:
        script=(tools_dir / "three_prime_seq" / "findSitesMM.py"),
        minlength=config["polyAT"]["minlength"],
        mindistance=config["polyAT"]["mindistance"],
        extension=config["polyAT"]["extend"],
        windowlength=config["polyAT"]["windowlength"],
        percbase=config["polyAT"]["percbase"],
    conda:
        CONDA_SHARED_ENV
    shell:
        "{params.script} -o {output} "
        "--tb {input.two_bit} "
        "--bed {input.bed} "
        "--minLength {params.minlength} "
        "--base {wildcards.base} "
        "--extend {params.extension} "
        "--minDistance {params.mindistance} "
        "--windowLength {params.windowlength} "
        "--percBase {params.percbase} "

# exclude second in pair, only 5' signal of R1
def bamcov_filter_opts(wc):
    s = ("--filterRNAstrand {direction} --samFlagExclude 128 "
         "--Offset 1 -bs 1 --skipNAs")
    return s.format(direction=wc.direction)

rule three_prime_seq_bam_cov:
    input:
        bam=aligner + "/{sample}.sorted.bam",
        bai=aligner + "/{sample}.sorted.bam.bai"
    output: 
        "three_prime_seq/raw/{sample}_direction-{direction}.bw"
    params:
        filterOpt=bamcov_filter_opts
    threads: 
        16
    conda:
        CONDA_SHARED_ENV
    shell: 
        "bamCoverage -b {input.bam} "
        "-p {threads} "
        "-o {output} "
        "{params.filterOpt} "


def filterbw_which_bed(wc):
    return ("three_prime_seq/polyA.bed" if wc.direction == "forward" 
            else "three_prime_seq/polyT.bed")

rule filterBW:
    input:
        bigwig="three_prime_seq/raw/{sample}_direction-{direction}.bw",
        bed=filterbw_which_bed
    output: 
        "three_prime_seq/filtered/{sample}_direction-{direction}.bw"
    params:
        script=(tools_dir / "three_prime_seq" / "filterBW.py")
    conda:
        CONDA_SHARED_ENV
    shell:
        "{params.script} {input} {output} "


# Associate signal with each gene (flank by some amount)
# Note how many bases couldn't be associated with any gene
rule geneAssociation:
    input:
        expand("three_prime_seq/filtered/{{sample}}_direction-{direction}.bw", direction=["forward", "reverse"])
    output: "three_prime_seq/{sample}_polyA_annotation.txt"
    threads: 16
    params:
        extension=config["geneAssociation"]["extend"], # 500
        gtf=genes_gtf,
        script=(tools_dir / "three_prime_seq" / "signal2gene.py")
    conda:
        CONDA_SHARED_ENV
    shell: 
        "{params.script} --extend {params.extension} "
        "--threads {threads} "
        "{input} {params.gtf} {output} "

# return list of replicates for each condition output by geneAssociation
# def find_replicates_cluster_pas(wc):
#     _samples = samples_from_condition(str(wc.condition))
#     return expand("three_prime_seq/{sample}_polyA_annotation.txt", sample=_samples)

# here we need to merge all replicates of geneAssociation -> see 
# /data/hilgers/group2/rezansoff/sakshiProj/2507_3seq/polyA_annotation_polysome2K/combined_cluster_andGrep_commands
# also sort by start position
rule preprocess_cluster_pas:
    input: 
        expand("three_prime_seq/{sample}_polyA_annotation.txt", sample=samples)
    output: 
        temp("three_prime_seq/tmp/geneassociation_merged_preprocessed.txt")
    shell:        
        "cat {input} | "
        "sed '/^[ ]*Chrom/ d' | " 
        "sort -k1,1 -k2,2n " 
        "> {output} "

rule clusterPAS: 
    # input is preprocessed output from geneAssoc
    input: 
        "three_prime_seq/tmp/geneassociation_merged_preprocessed.txt"
    output: 
        temp("three_prime_seq/tmp/clusterPAS_tmpdb.txt")
    conda:
        CONDA_SHARED_ENV
    params: 
        script=(tools_dir / "three_prime_seq" / "clusterPAS.py"),
        windowsize=config["clusterPAS"]["window"], # 15
    shell:
        "python {params.script} --windowSize {params.windowsize} {input} {output}"


# awk command: remove entries with multiple genes in 4th column (must be unambiguous)
# python script: add "_1", "_2", to each cluster label (4th column) to make each 
# unique for each genomic position
# TODO possibly make filtering of CDS/exons optional
# TODO possibly keep only those that intersect an annotated 3' UTR
rule postprocess_cluster_pas:
    input: 
        "three_prime_seq/tmp/clusterPAS_tmpdb.txt"
    output: 
        "three_prime_seq/db/clusterPAS_db.txt"
    conda:
        CONDA_SHARED_ENV
    params:
        dedup_script=(tools_dir / "three_prime_seq" / "dedup_clusterPAS.py")
    shell:
        """
        grep -v "5'UTR" {input} | \
        grep -v "CDS" | \
        grep -v "exon" | \
        awk '!($4 ~ /[,]/)' | \
        python {params.dedup_script} > {output}
        """


# return clusterPAS db for corresponding condition given sample
# def find_cluster_pas_db(wc):
#     condition = condition_from_sample(str(wc.sample))
#     return ("three_prime_seq/db/condition-{condition}_clusterPAS_db.txt"
#             .format(condition=condition))


# previously done by hand 
# e.g. /data/hilgers/group2/rezansoff/sakshiProj/2507_3seq/countReadEnds_commands_real.txt
# each line of output of clusterPAS must be unique'd 
rule count_read_ends:
    input: 
        bws=expand("three_prime_seq/filtered/{{sample}}_direction-{direction}.bw", 
                   direction=["forward", "reverse"]),
        bed="three_prime_seq/db/clusterPAS_db.txt", # this is the "database" of PAS sites
    output:
        counts="three_prime_seq/{sample}_uniqcounts.txt"
    wildcard_constraints:
        sample="[^\/]+" # no /
    conda:
        CONDA_SHARED_ENV
    params:
        script=(tools_dir / "three_prime_seq" / "countReadEnds.py")
    shell:
        "python {params.script} {input} {output}"

# create counts.tsv of all samples for each cluster
rule merge_read_ends:
    input:
        expand("three_prime_seq/{sample}_uniqcounts.txt", sample=samples)
    output: 
        "three_prime_seq/counts.tsv"
    conda:
        CONDA_SHARED_ENV
    params:
        script=(tools_dir / "three_prime_seq" / "mergeReadEnds.py"),
        samples=samples
    shell:
        "python {params.script} -i {input} -s {params.samples} -o {output}"


rule cmatrix_raw:
    input: 
        expand("three_prime_seq/raw/{sample}_direction-{{direction}}.bw", sample=samples)
    output: 
        temp("three_prime_seq/cmatrix_raw_direction-{direction}.mat.gz")
    threads: 
        16
    conda:
        CONDA_SHARED_ENV
    params:
        upstream=config["cmatrix_raw"]["upstream"], #500
        downstream=config["cmatrix_raw"]["downstream"], #500
        labels=samples,
        bed=genes_bed,
    shell: 
        "computeMatrix scale-regions "
        "-S {input} "
        "-R {params.bed} "
        "--metagene "
        "--unscaled3prime {params.upstream} "
        "-a {params.downstream} "
        "-p {threads} "
        "--samplesLabel {params.labels} "
        "--skipZeros "
        "-bs 5 "
        "-o {output} "


rule cmatrix_filtered:
    input: 
        "three_prime_seq/cmatrix_raw_direction-{direction}.mat.gz"
    output: 
        temp("three_prime_seq/cmatrix_filtered_direction-{direction}.mat.gz")
    threads: 8
    params:
        strand=lambda wc: "+" if wc.direction == "forward" else "-"
    conda: 
        CONDA_SHARED_ENV
    shell: 
        "computeMatrixOperations filterStrand "
        "-m {input} "
        "-o {output} "
        "--strand '{params.strand}' "


rule merge_matrix:
    input:
        expand("three_prime_seq/cmatrix_filtered_direction-{direction}.mat.gz",
               direction=["forward", "reverse"])
    output: 
        "three_prime_seq/combined_polyA.mat.gz"
    threads: 8
    conda: 
        CONDA_SHARED_ENV
    shell: 
        "computeMatrixOperations rbind "
        "-m {input} "
        "-o {output} "

rule heatmap:
    input: "three_prime_seq/combined_polyA.mat.gz"
    output: "three_prime_seq/combined_polyA.png"
    conda: 
        CONDA_SHARED_ENV
    shell: 
        "plotProfile -m {input} -o {output}"

if sampleSheet:
    rule DESeq2:
        input:
            gene_counts="featureCounts/counts.tsv",
            cluster_counts="three_prime_seq/counts.tsv",
            sampleSheet=sampleSheet,
            symbol_file = "Annotation/genes.filtered.symbol"
        output:
            "{}/DESeq2.session_info.txt".format(get_outdir("DESeq2",sampleSheet))
        benchmark:
            "{}/.benchmark/DESeq2.featureCounts.benchmark".format(get_outdir("DESeq2",sampleSheet))
        params:
            outdir = os.path.join(outdir, get_outdir("DESeq2", sampleSheet)),
            fdr = 0.05,
        conda: CONDA_RNASEQ_ENV
        script: "../rscripts/threeprimeseq-DESeq2.R"
