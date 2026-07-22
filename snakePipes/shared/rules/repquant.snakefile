import pandas as pd

df=pd.read_csv(sampleSheet,sep="\t")
if isMultipleComparison:
    combined_col="condition_group"
    df[combined_col] = df["condition"] + "_" + df["group"]
    groupby_col=combined_col
else:
    groupby_col="condition"
merge_dict=df.groupby(groupby_col)["name"].apply(list).to_dict()


rule merge_peaks:
    input:
        peaks=lambda wildcards: expand("SEACR/{chip_sample}.filtered.stringent.bed", chip_sample=merge_dict[wildcards.merge_group]),
        sampleSheet = sampleSheet
    output:
        merged_peaks = "SEACR_merged_peaks/{merge_group}_peak_intersect_f0.5.bed"
    params:
        b = lambda wildcards,input: " ".join(f"-b {sfile}" for sfile in input.peaks[1:])
    conda: CONDA_SEACR_ENV
    shell: """
            bedtools intersect -a {input.peaks[0]} {params.b} -f 0.5 > {output}
           """

rule makeRMSKGTF:
    input: rmsk_file
    output: "Annotation/rmsk.gtf"
    run:
        f = open(input[0])
        of = open(output[0], "w")
        found = dict()

        for line in f:
            if line.startswith("#"):
                continue
            bin, swScore, milliDiv, milliDel, milliIns,  genoName, genoStart, genoEnd, genoLeft, strand, repName, repClass, repFamily, repStart, repEnd, repLeft, id = line.strip().split("\t")

            genoStart_fixed = int(genoStart) + 1 #standard gtf uses 1-based coordinates

            if repName not in found:
                found[repName] = 0
                tid = repName
            else:
                found[repName] += 1
                tid = "{}_dup{}".format(repName, found[repName])

            meta = {"genoName": genoName,
                    "genoStart": genoStart_fixed,
                    "genoEnd": genoEnd,
                    "strand": strand,
                    "repName": repName,
                    "repClass": repClass,
                    "repFamily": repFamily,
                    "tid": tid}

            of.write("{genoName}\trmsk\texon\t{genoStart}\t{genoEnd}\t.\t{strand}\t.\tgene_id \"{repName}\"; transcript_id \"{tid}\"; family_id \"{repFamily}\"; class_id \"{repClass}\";\n".format(**meta))
        f.close()
        of.close()


rule intersect_peaks_rmsk:
    input:
        rmsk_gtf="Annotation/rmsk.gtf",
        merged_peaks = "SEACR_merged_peaks/{merge_group}_peak_intersect_f0.5.bed"
    output:
        rmsk_annotated_peaks = "SEACR_annotated_peaks/{merge_group}_rmsk.bed"
    conda: CONDA_SEACR_ENV
    shell: """
           bedtools intersect -wa -u -a <( bedtools sort -i {input.rmsk_gtf} ) -b <( bedtools sort -i {input.merged_peaks} ) > {output.rmsk_annotated_peaks}
           """        

rule randomize_peaks:
    input:
        merged_peaks = "SEACR_merged_peaks/{merge_group}_peak_intersect_f0.5.bed"
    output:
        randomized_peaks = "SEACR_randomized_peaks/{merge_group}_randomized_peaks.bed"
    params:
        script = os.path.join(maindir, "shared", "rscripts","repquant_rand_ranges.R"),
        outdir = "SEACR_randomized_peaks",
        merged_peaks = lambda wildcards,input: os.path.join(outdir,input.merged_peaks)
    conda: CONDA_REPQUANT_ENV
    script: "{params.script}"


rule intersect_randomized_peaks_rmsk:
    input:
        rmsk_gtf="Annotation/rmsk.gtf",
        randomized_peaks = "SEACR_randomized_peaks/{merge_group}_randomized_peaks.bed"
    output:
        rmsk_annotated_peaks = "SEACR_randomized_peaks/{merge_group}_rmsk.bed"
    conda: CONDA_SEACR_ENV
    shell: """
           bedtools intersect -wa -u -a <( bedtools sort -i {input.rmsk_gtf} ) -b <( bedtools sort -i {input.randomized_peaks} ) > {output.rmsk_annotated_peaks}
           """


rule generate_repquant_report:
    input:
        rmsk_annotated_peaks = expand("SEACR_annotated_peaks/{merge_group}_rmsk.bed",merge_group=merge_dict.keys()),
        rmsk_annotated_randomized_peaks = expand("SEACR_randomized_peaks/{merge_group}_rmsk.bed",merge_group=merge_dict.keys()),
        rmsk_gtf="Annotation/rmsk.gtf"
    output:
        report_html = "RepQuant/report.html"
    params:
        script = os.path.join(maindir, "shared", "rscripts","repquant_report.Rmd"),
        outdir = "RepQuant",
        rmsk_annotated_peaks = lambda wildcards,input: [os.path.join(outdir,x) for x in input.rmsk_annotated_peaks],
        rmsk_annotated_randomized_peaks = lambda wildcards,input: [os.path.join(outdir,x) for x in input.rmsk_annotated_randomized_peaks],
        rmsk_gtf=lambda wildcards,input: os.path.join(outdir,input.rmsk_gtf)
    conda: CONDA_REPQUANT_ENV
    script: "{params.script}"
