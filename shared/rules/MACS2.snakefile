def generate_bam_pairs_for_MACS2(wildcards):
    chip = wildcards.sample
    inpt = None # corresponding sample name
    for g in groups:
        if chip == g[0]:
            inpt = g[1]
            break
    ##print(chip, inpt)
    return( expand("MACS2/{sample}.noDuplicates.bam", sample=(chip,inpt)) )


rule noDuplicates_for_MACS2:
    input:  "Bowtie2/{sample}.bam"
    output: "MACS2/{sample}.noDuplicates.bam"
    benchmark:  "MACS2/.benchmark/noDuplicates_for_MACS2.{sample}.benchmark"
    shell:  samtools_path+" view -hb -F1024 {input} > {output}"


if paired:
    rule MACS2:
        input:
            bams = generate_bam_pairs_for_MACS2,
            mean_file = "InsertSizeMetrics/{sample}.mean.txt"
        output: "MACS2/{sample}_peaks.xls"
        params:
            description = "{sample}",
            fragment_length = fragment_length,
            genome_size = genome_size
        log:    "MACS2/log/{sample}.log"
        benchmark:  "MACS2/.benchmark/MACS2.{sample}.benchmark"
        run:
            try:
                fragment_length = int(get_from_file("mean", input["mean_file"]))
            except:
                fragment_length = default_fragment_length
            input_str = " -c ".join(input["bams"])
            shell(
                "source "+macs2_activate+" && "
                ""+macs2_path+" callpeak --nomodel --extsize {fragment_length} "
                "-n {params.description} --gsize {params.genome_size} "
                "-t {input_str} "
                "--outdir MACS2 "
                #"&& rm {input} "   # not possible to remove files here when one input is used for multiple ChIPs!!!
                "2>&1 | tee {log}")
else:
    rule MACS2:
        input:  generate_bam_pairs_for_MACS2
        output: "MACS2/{sample}_peaks.xls"
        params:
            description = "{sample}",
            fragment_length = fragment_length,
            genome_size = genome_size
        log:    "MACS2/log/{sample}.log"
        benchmark:  "MACS2/.benchmark/MACS2.{sample}.benchmark"
        run:
            input_str = " -c ".join(input)
            shell(
                "source "+macs2_activate+" && "
                ""+macs2_path+" callpeak --nomodel --extsize {params.fragment_length} "
                "-n {params.description} --gsize {params.genome_size} "
                "-t {input_str} "
                "--outdir MACS2 "
                #"&& rm {input} "   # not possible to remove files here when one input is used for multiple ChIPs!!!
                "2>&1 | tee {log}")
