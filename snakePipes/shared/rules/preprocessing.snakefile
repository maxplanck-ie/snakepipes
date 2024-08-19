# This requires sampleDict, which is a dictionary defined in the preprocessing Snakefile
import os

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
            run:
                if not os.path.exists(os.path.join(outdir,output.r1)):
                    os.symlink(os.path.join(outdir,input.r1),os.path.join(outdir,output.r1))
                if not os.path.exists(os.path.join(outdir,output.r2)):
                    os.symlink(os.path.join(outdir,input.r2),os.path.join(outdir,output.r2))
    else:
        rule mergeFastq:
            input:
                r1=indir + "/{sample}" + reads[0] + ext
            output:
                r1="mergedFASTQ/{sample}" + reads[0] + ext
            run:
                if not os.path.exists(os.path.join(outdir,output.r1)):
                    os.symlink(os.path.join(outdir,input.r1),os.path.join(outdir,output.r1))

if optDedupDist > 0:
    if pairedEnd:
        rule clumpify:
            input:
                r1="mergedFASTQ/{sample}" + reads[0] + ext,
                r2="mergedFASTQ/{sample}" + reads[1] + ext
            output:
                r1="deduplicatedFASTQ/{sample}" + reads[0] + ext,
                r2="deduplicatedFASTQ/{sample}" + reads[1] + ext,
                metrics="deduplicatedFASTQ/{sample}.metrics",
                tempOut=temp("deduplicatedFASTQ/{sample}_temp.fq.gz")
            params:
                R1=reads[0],
                R2=reads[1],
                extension=ext,
                mem=clumpifyMemory,
                optDedupDist=optDedupDist,
                clumpifyOptions=clumpifyOptions
            benchmark: "deduplicatedFASTQ/.benchmarks/{sample}"
            threads: lambda wildcards: 20 if 20<max_thread else max_thread
            conda: CONDA_PREPROCESSING_ENV
            shell: """
                clumpify.sh -Xmx{params.mem} \
                            {params.clumpifyOptions} \
                            in={input.r1} \
                            in2={input.r2} \
                            out={output.tempOut} \
                            dupedist={params.optDedupDist} \
                            threads={threads}

                splitFastq --pigzThreads 4 --R1 {params.R1} --R2 {params.R2} --extension {params.extension} \
                           {output.tempOut} \
                           deduplicatedFASTQ/{wildcards.sample} > {output.metrics}
                """
    else:
        rule clumpify:
            input:
                r1="mergedFASTQ/{sample}" + reads[0] + ext,
            output:
                r1="deduplicatedFASTQ/{sample}" + reads[0] + ext,
                metrics="deduplicatedFASTQ/{sample}.metrics",
                tempOut=temp("deduplicatedFASTQ/{sample}_temp.fq.gz")
            params:
                R1=reads[0],
                extension=ext,
                mem=clumpifyMemory,
                optDedupDist=optDedupDist,
                clumpifyOptions=clumpifyOptions
            benchmark: "deduplicatedFASTQ/.benchmarks/{sample}"
            threads: lambda wildcards: 20 if 20<max_thread else max_thread
            conda: CONDA_PREPROCESSING_ENV
            shell: """
                clumpify.sh -Xmx{params.mem} \
                            {params.clumpifyOptions} \
                            in={input.r1} \
                            out={output.tempOut} \
                            dupedist={params.optDedupDist} \
                            threads={threads}

                splitFastq --SE --pigzThreads 4 --R1 {params.R1} --extension {params.extension} \
                           {output.tempOut} \
                           deduplicatedFASTQ/{wildcards.sample} > {output.metrics}
                """
else:
    if pairedEnd:
        rule clumpify:
            input:
                r1="mergedFASTQ/{sample}" + reads[0] + ext,
                r2="mergedFASTQ/{sample}" + reads[1] + ext
            output:
                r1="deduplicatedFASTQ/{sample}" + reads[0] + ext,
                r2="deduplicatedFASTQ/{sample}" + reads[1] + ext
            shell: """
                ln -s ../{input.r1} {output.r1};
                ln -s ../{input.r2} {output.r2}
            """
    else:
        rule clumpify:
            input:
                r1="mergedFASTQ/{sample}" + reads[0] + ext
            output:
                r1="deduplicatedFASTQ/{sample}" + reads[0] + ext
            shell: """
                ln -s ../{input.r1} {output.r1}
            """

rule splitFastq2YAML:
    input:
        expand("deduplicatedFASTQ/{sample}.metrics", sample=samples)
    output:
        "deduplicatedFASTQ/optical_dedup_mqc.json"
    run:
        import json

        d = {
             "id": "custom_barplot",
             "section_name": "Optical Duplicates",
             "description": "The percentage of optical duplicates found in the original fastq files.",
             "plot_type": "bargraph",
             "pconfig": {
                         "id": "custom_bargraph",
                         "ylab": "Percentage optical duplicates",
                        },
             "data": {}
            }

        for fname in input:
            # name things according to the sample
            sample = fname.split("/")[1][:-8]
            f = open(fname)
            cols = f.read().strip().split("\t")
            cols = [int(x) for x in cols]
            if cols[1] == 0:
                cols[1] = 1
            d["data"][sample] = {"Optical Duplicates": cols[0],
                                 "Non-Duplicates": cols[1]}
            f.close()

        o = open(output[0], "w")
        json.dump(d, o)
        o.close()
