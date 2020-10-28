## function to get the name of the samplesheet and extend the name of the folder. (Same logic as for DESeq2).
def get_outdir(folder_name,sampleSheet):
    sample_name = os.path.splitext(os.path.basename(str(sampleSheet)))[0]
    return("{}_{}".format(folder_name, sample_name))
## reWrap libType to string to use as flag for rmats rather than int.
def wrap_libType(libType):
    dic_libType = {0:"fr-unstranded",1:"fr-firststrand",2:"fr-secondstrand"}
    return dic_libType[libType]

## requires the checkpoint rule defined in DESeq2.multipleComp.snakefile
#rMatsConds = cf.sampleSheetGroups(sampleSheet)

def generate_b1_b2(sampleSheet,which_b):
    rMatsConds = cf.sampleSheetGroups(sampleSheet)
    if which_b == "b1":
        return ",".join(["filtered_bam/" + s for s in [s + ".filtered.bam" for s in rMatsConds[list(rMatsConds)[0]]]])
    else:
        return ",".join(["filtered_bam/" + s for s in [s + ".filtered.bam" for s in rMatsConds[list(rMatsConds)[1]]]])

def get_s1(sampleSheet):
    rMatsConds = cf.sampleSheetGroups(sampleSheet)
    return ["filtered_bam/" + s for s in [s + ".filtered.bam" for s in rMatsConds[list(rMatsConds)[0]]]][0]

rule createInputcsv:
    input:
        bams = expand("filtered_bam/{sample}.filtered.bam.bai", sample=samples),
        sampleSheet = "splitSampleSheets/" + os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv"
    output:
        b1out = "rMats_{}/b1.csv".format(os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}"),
        b2out = "rMats_{}/b2.csv".format(os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}")
    params:
        b1 = lambda wildcards,input: generate_b1_b2(os.path.join(outdir,str(input.sampleSheet)),"b1"),
        b2 = lambda wildcards,input: generate_b1_b2(os.path.join(outdir,str(input.sampleSheet)),"b2")
    shell: """
        echo '{params.b1}' > {output.b1out}
        echo '{params.b2}' > {output.b2out}
"""

rule rMats:
    input:
        b1 = "rMats_{}/b1.csv".format(os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}"),
        b2 = "rMats_{}/b2.csv".format(os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}"),
        sampleSheet = "splitSampleSheets/" + os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}.tsv"
    output:
        "rMats_{}/RI.MATS.JCEC.txt".format(os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}")
    params:
        s1 = lambda wildcards,input: get_s1(str(input.sampleSheet)),
        readLen = "rMats_{}/readlength.txt".format(os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}"),
        gtf = genes_gtf,
        od = "rMats_{}".format(os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}"),
        end = "paired" if pairedEnd else "single",
        libType = wrap_libType(libraryType),
        tempDir = tempDir,
    log: "rMats_{}/rMats.log".format(os.path.splitext(os.path.basename(str(sampleSheet)))[0]+".{compGroup}")
    threads: 4
    conda: CONDA_RMATS_ENV
    shell:"""
        TMPDIR={params.tempDir}
        MYTEMP=$(mktemp -d ${{TMPDIR:-/tmp}}/snakepipes.XXXXXXXXXX);
        set +o pipefail;
        readLen=$(samtools view {params.s1} | awk \'{{print length($10)}}\' | head -10000 | awk \'{{ sum += $1 }} END {{ if (NR > 0) print sum / NR }}\')
        rmats.py --gtf {params.gtf} --b1 {input.b1} --b2 {input.b2} --od {params.od} --tmp $MYTEMP -t {params.end} --libType {params.libType} --readLength $readLen --variable-read-length --nthread {threads} --tstat {threads} 2> {log};
        rm -rf $MYTEMP
"""
