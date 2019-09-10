# This requires sampleDict, which is a dictionary defined in the preprocessing Snakefile

initialIndir = indir  # Due to the lambda functions, this ends up getting reset over time
if sampleSheet:
    if pairedEnd:
        rule mergeFastq:
            input:
                r1=lambda wildcards: expand(initialIndir + "/{sample}", sample=sampleDict[wildcards.sample][0]),
                r2=lambda wildcards: expand(initialIndir + "/{sample}", sample=sampleDict[wildcards.sample][1])
            output:
                r1="mergedFASTQ/{sample}" + reads[0] + ext,
                r2="mergedFASTQ/{sample}" + reads[1] + ext
            shell: """
                cat {input.r1} > {output.r1}
                cat {input.r2} > {output.r2}
                """
    else:
        print("SE")
        rule mergeFastq:
            input:
                r1=lambda wildcards: expand(initialIndir + "/{sample}", sample=sampleDict[wildcards.sample][0])
            output:
                r1="mergedFASTQ/{sample}" + reads[0] + ext
            shell: """
                cat {input.r1} > {output.r1}
                """
else:
    if pairedEnd:
        rule mergeFastq:
            input:
                r1=indir + "/{sample}" + reads[0] + ext,
                r2=indir + "/{sample}" + reads[1] + ext
            output:
                r1="mergedFASTQ/{sample}" + reads[0] + ext,
                r2="mergedFASTQ/{sample}" + reads[1] + ext
            shell: """
                ln -sr {input.r1} {output.r1}
                ln -sr {input.r2} {output.r2}
                """
    else:
        rule mergeFastq:
            input:
                r1=indir + "/{sample}" + reads[0] + ext
            output:
                r1="mergedFASTQ/{sample}" + reads[0] + ext
            shell: """
                ln -sr {input.r1} {output.r1}
                """

if optDedupDist > 0:
    if pairedEnd:
        rule clumpify:
            input:
                r1="mergedFASTQ/{sample}" + reads[0] + ext,
                r2="mergedFASTQ/{sample}" + reads[1] + ext
            output:
                r1="deduplicatedFASTQ/{sample}" + reads[0] + ext,
                r2="deduplicatedFASTQ/{sample}" + reads[1] + ext,
                metrics="deduplicatedFASTQ/{sample}.metrics"
            params:
                R1=reads[0],
                R2=reads[1],
                extension=ext,
                mem=clumpifyMemory,
                optDedupDist=optDedupDist,
                clumpifyOptions=clumpifyOptions
            log:
                stdout="deduplicatedFASTQ/logs/{sample}.stdout",
                stderr="deduplicatedFASTQ/logs/{sample}.stderr"
            benchmark: "deduplicatedFASTQ/.benchmarks/{sample}"
            threads: 20
            conda: CONDA_PREPROCESSING_ENV
            shell: """
                clumpify.sh -Xmx{params.mem} \
                            {params.clumpifyOptions} \
                            in={input.r1} \
                             in2={input.r2} \
                            out=deduplicatedFASTQ/{wildcards.sample}_temp.fq.gz \
                            dupedist={params.optDedupDist} \
                            threads={threads} > {log.stdout} 2> {log.stderr}

                splitFastq --pigzThreads 4 --R1 {params.R1} --R2 {params.R2} --extension {params.extension} \
                           deduplicatedFASTQ/{wildcards.sample}_temp.fq.gz \
                           deduplicatedFASTQ/{wildcards.sample} > {output.metrics} 2>> {log.stderr}

                rm deduplicatedFASTQ/{wildcards.sample}_temp.fq.gz
                """
    else:
        rule clumpify:
            input:
                r1="mergedFASTQ/{sample}" + reads[0] + ext,
            output:
                r1="deduplicatedFASTQ/{sample}" + reads[0] + ext,
                metrics="deduplicatedFASTQ/{sample}.metrics"
            params:
                R1=reads[0],
                extension=ext,
                mem=clumpifyMemory,
                optDedupDist=optDedupDist,
                clumpifyOptions=clumpifyOptions
            log:
                stdout="deduplicatedFASTQ/logs/{sample}.stdout",
                stderr="deduplicatedFASTQ/logs/{sample}.stderr"
            benchmark: "deduplicatedFASTQ/.benchmarks/{sample}"
            threads: 20
            conda: CONDA_PREPROCESSING_ENV
            shell: """
                clumpify.sh -Xmx{params.mem} \
                            {params.clumpifyOptions} \
                            in={input.r1} \
                            out=deduplicatedFASTQ/{wildcards.sample}_temp.fq.gz \
                            dupedist={params.optDedupDist} \
                            threads={threads} > {log.stdout} 2> {log.stderr}

                splitFastq --SE --pigzThreads 4 --R1 {params.R1} --extension {params.extension} \
                           deduplicatedFASTQ/{wildcards.sample}_temp.fq.gz \
                           deduplicatedFASTQ/{wildcards.sample} > {output.metrics} 2>> {log.stderr}

                rm deduplicatedFASTQ/{wildcards.sample}_temp.fq.gz
                """
else:
    if pairedEnd:
        rule clumpify:
            input:
                r1="mergedFASTQ/{sample}" + reads[0] + ext,
                r2="mergedFASTQ/{sample}" + reads[1] + ext
            output:
                r1="deduplicatedFASTQ/{sample}" + reads[0] + ext,
                r2="deduplicatedFASTQ/{sample}" + reads[1] + ext,
            shell: """
                ln -s -r {input.r1} {output.r1}
                ln -s -r {input.r2} {output.r2}
                """
    else:
        rule clumpify:
            input:
                r1="mergedFASTQ/{sample}" + reads[0] + ext
            output:
                r1="deduplicatedFASTQ/{sample}" + reads[0] + ext
            shell: """
                ln -s -r {input.r1} {output.r1}
                """
