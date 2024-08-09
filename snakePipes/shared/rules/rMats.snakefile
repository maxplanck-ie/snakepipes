## function to get the name of the samplesheet and extend the name of the folder. (Same logic as for DESeq2).
def get_outdir(folder_name,sampleSheet):
    sample_name = os.path.splitext(os.path.basename(str(sampleSheet)))[0]
    return("{}_{}".format(folder_name, sample_name))
## reWrap libType to string to use as flag for rmats rather than int.
def wrap_libType(libType):
    dic_libType = {0:"fr-unstranded",1:"fr-firststrand",2:"fr-secondstrand"}
    return dic_libType[libType]

rule createInputcsv:
    input:
        expand("filtered_bam/{sample}.filtered.bam.bai", sample=samples)
    output:
        b1out = "{}/b1.csv".format(get_outdir("rMats", sampleSheet)),
        b2out = "{}/b2.csv".format(get_outdir("rMats", sampleSheet))
    params:
        b1 = ",".join(["filtered_bam/" + s for s in [s + ".filtered.bam" for s in rMatsConds[list(rMatsConds)[0]]]]),
        b2 = ",".join(["filtered_bam/" + s for s in [s + ".filtered.bam" for s in rMatsConds[list(rMatsConds)[1]]]])
    shell: """
        echo '{params.b1}' > {output.b1out}
        echo '{params.b2}' > {output.b2out}
"""
rule rMats:
    input:
        b1 = "{}/b1.csv".format(get_outdir("rMats", sampleSheet)),
        b2 = "{}/b2.csv".format(get_outdir("rMats", sampleSheet))
    output:
        "{}/RI.MATS.JCEC.txt".format(get_outdir("rMats", sampleSheet))
    params:
        s1 = ["filtered_bam/" + s for s in [s + ".filtered.bam" for s in rMatsConds[list(rMatsConds)[0]]]][0],
        readLen = "{}/readlength.txt".format(get_outdir("rMats", sampleSheet)),
        gtf = genes_gtf,
        od = "{}".format(get_outdir("rMats", sampleSheet)),
        end = "paired" if pairedEnd else "single",
        libType = wrap_libType(libraryType),
        tempDir = tempDir,
    threads: 1
    conda: CONDA_RNASEQ_ENV
    shell:"""
        set +o pipefail;
        readLen=$(samtools view {params.s1} | awk \'{{print length($10)}}\' | head -10000 | awk \'{{ sum += $1 }} END {{ if (NR > 0) print int(sum / NR) }}\')
        rmats.py --gtf {params.gtf} --b1 {input.b1} --b2 {input.b2} --od {params.od} --tmp {params.tempDir} -t {params.end} --libType {params.libType} --readLength $readLen --variable-read-length --nthread {threads} --tstat {threads}
"""
