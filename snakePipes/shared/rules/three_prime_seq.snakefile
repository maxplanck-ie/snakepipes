# adapted from Andrew Rezansoff's 3' seq pipeline tools
# /data/hilgers/group/rezansoff/3seq_pipeline_tools/

# TODO: 
# 1. pass three_prime_seq fastp arguments to fastp
# 2. specify all output in mRNA-seq snakefile -> how to pass --three-prime-seq?
# 3. check why cmatrix_filtered is commented out

tools_dir = Path(maindir) / "shared" / "tools"

# trimming done using STAR and fastp

rule all_three_prime_samples:
    input: 

rule clusterPAS:
    input: 
        "three_prime_seq/SampleAll.txt"
    conda:
        CONDA_SHARED_ENV
    params: 
        script=(tools_dir / "three_prime_seq" / "clusterPAS.py"),
        windowsize=config["clusterPASWindowSize"], # 10

rule polyAT: 
    input: 
        twobit=genome_2bit,
        bed=genes_bed
    output: 
        "three_prime_seq/poly{base}.bed"
    params:
        script=(tools_dir / "three_prime_seq" / "findSitesMM.py"),
        minlength=config["polyAT"]["minlength"],
        mindistance=config["polyAT"]["mindistance"],
        extend=config["polyAT"]["extend"],
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
        "--extend {params.extend} "
        "--minDistance {params.minDistance} "
        "--windowLength {params.windowLength} "
        "--percBase {params.percBase} "

def bamcov_filter_opts(wc):
    s = ("--filterRNAstrand {direction} --samFlagExclude 128 
         "--Offset 1 -bs 1 --skipNAs")
    return s.format(direction=wc.direction)

rule three_prime_seq_bam_cov:
    input:
        bam=aligner + "/{sample}.sorted.bam",
        bai=aligner + "/{sample}.sorted.bam.bai"
    output: 
        "three_prime_seq/{sample}_{direction}.bw"
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
    return ("three_prime_seq/polyA.bed" if wildcards.direction == "forward" 
            else "three_prime_seq/polyT.bed")

rule filterBW:
    input:
        bigwig="three_prime_seq/{sample}_{direction}.bw",
        bed=filterbw_which_bed
    output: 
        "three_prime_seq/{sample}_{direction}.filtered.bw"
    log:
        stdout="three_prime_seq/logs/{sample}_{direction}.stdout,
        stderr="three_prime_seq/logs/{sample}_{direction}.stderr"
    params:
        script=(tools_dir / "three_prime_seq" / "filterBW.py")
    conda:
        CONDA_SHARED_ENV
    shell:
        "{params.script} {input} {output} "
        "> {log.stdout} "
        "2> {log.stderr} "

# previously done by hand 
# e.g. /data/hilgers/group2/rezansoff/sakshiProj/2507_3seq/countReadEnds_commands_real.txt
rule count_read_ends:
    input: 
        expand("three_prime_seq/{{sample}}_{direction}.filtered.bw", 
               direction=["forward", "reverse"])
    output:
        ids="three_prime_seq/{sample}_uniqIDs.txt",
        counts="three_prime_seq/{sample}_uniqcounts.txt"
    params:
        script=(tools_dir / "three_prime_seq / "countReadEnds.py")
    shell:
        "{params.script} {input} {output}"

# Associate signal with each gene (flank by some amount)
# Note how many bases couldn't be associated with any gene
rule geneAssociation:
    input:
        expand("three_prime_seq/{{sample}}_{direction}.filtered.bw", 
               direction=["forward", "reverse"])
    output: "three_prime_seq/{sample}_polyA_annotation.txt"
    threads: 16
    params:
        extend=config["geneAssociation"]["extend"], # 2000
        gtf=genes_gtf,
        script=(tools_dir / "three_prime_seq" / "signal2gene.py")
    conda:
        CONDA_SHARED_ENV
    shell: 
        "{params.script} --extend {params.extend} "
        "--threads {threads} "
        "{input} {params.gtf} {output} "


## Need a better regex filter ".+(?<!\.filtered)"
rule cmatrix_raw:
    input: 
        expand("three_prime_seq/{sample}_{{direction}}.bw", sample = samples)
    output: 
        temp("three_prime_seq/cmatrix_raw_direction-{direction}.mat.gz")
    threads: 
        16
    conda:
        CONDA_SHARED_ENV
    params:
        upstream=config["geneAssociation"]["upstream"], #500
        downstream=config["geneAssociation"]["downstream"], #500
        labels=samples,
        bed=genes_bed
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
        temp("three_prime_seq/cmatrix_raw_direction-{direction}.filtered.mat.gz")
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
        expand("three_prime_seq/cmatrix_raw_direction-{direction}.filtered.mat.gz",
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