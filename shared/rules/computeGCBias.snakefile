if paired:
    rule computeGCBias:
        input:
            bam = "Bowtie2/{sample}.bam",
            index = "Bowtie2/{sample}.bam.bai",
            genome_2bit = genome_2bit,
            fragment_length_file = "InsertSizeMetrics/{sample}.mean.txt"
        output:
            plot = "GCBias/{sample}.GCBias.png",
            tsv = "GCBias/{sample}.GCBias.tsv"
        params:
            genome_size = genome_size,
            sample_size = 1000000    # must be a float!
        log:    "GCBias/log/{sample}.log"
        benchmark:  "GCBias/.benchmark/deepTools_computeGCBias.{sample}.benchmark"
        threads: 20
        shell:  os.path.join(deepTools2_dir,"computeGCBias") + " -p {threads} --sampleSize {params.sample_size} "
                    "-b {input.bam} --genome {input.genome_2bit} "
                    "--effectiveGenomeSize {params.genome_size} "
                    "--biasPlot {output.plot} --GCbiasFrequenciesFile {output.tsv} "
                    "2>&1 | tee {log} "
else:
    rule computeGCBias:
        input:
            bam = "Bowtie2/{sample}.bam",
            index = "Bowtie2/{sample}.bam.bai",
            genome_2bit = genome_2bit
        output:
            plot = "GCBias/{sample}.GCBias.png",
            tsv = "GCBias/{sample}.GCBias.tsv"
        params:
            fragment_length = fragment_length,
            genome_size = genome_size,
            sample_size = 1000000    # must be a float!
        log:    "GCBias/log/{sample}.log"
        benchmark:  "GCBias/.benchmark/deepTools_computeGCBias.{sample}.benchmark"
        threads: 20
        shell:  os.path.join(deepTools2_dir,"computeGCBias") + " -p {threads} --sampleSize {params.sample_size} "
                    "-b {input.bam} --genome {input.genome_2bit} "
                    "--effectiveGenomeSize {params.genome_size} --fragmentLength {params.fragment_length} "
                    "--biasPlot {output.plot} --GCbiasFrequenciesFile {output.tsv} "
                    "2>&1 | tee {log} "
