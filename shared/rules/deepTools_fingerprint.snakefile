def generate_paired_bam_from_aligned_files(wildcards):
    chip = wildcards.sample
    inpt = None # corresponding input name
    for g in groups:
        if chip == g[0]:
            inpt = g[1]
            break
    ##print(chip, inpt)
    return( expand("Bowtie2/{sample}.bam", sample=(chip,inpt)) )


def generate_paired_bai_from_aligned_files(wildcards):
    chip = wildcards.sample
    inpt = None # corresponding input name
    for g in groups:
        if chip == g[0]:
            inpt = g[1]
            break
    ##print(chip, inpt)
    return( expand("Bowtie2/{sample}.bam.bai", sample=(chip,inpt)) )


rule plotFingerprint:
    input:
        bam = generate_paired_bam_from_aligned_files,
        bai = generate_paired_bai_from_aligned_files
    output: "plotFingerprint/{sample}.plotFingerprint.png"
    params:
        binsize = 25,
        fragment_length = default_fragment_length     # Relevant for single-end paired_read_ext only!
    log:    "plotFingerprint/logs/{sample}.log"
    benchmark:  "plotFingerprint/.benchmark/deepTools_plotFingerprint.{sample}.benchmark"
    threads : 30
    run:
        input = " ".join(input["bam"])
        ##print("joined input:", input)
        shell(  os.path.join(deepTools2_dir,"plotFingerprint") + " -p {threads} "
                "-b {input} "
                "--binSize={params.binsize} "
                "--extendReads {params.fragment_length} "
                "--plotFile {output} "
                "2>&1 | tee {log} " )
