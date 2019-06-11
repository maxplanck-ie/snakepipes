import os
import re
from operator import is_not
import tempfile

###symlink bams if this is the starting point
if fromBAM:
    rule link_bam:
        input:
            indir + "/{sample}" + bamExt
        output:
            "bwameth/{sample}.PCRrm.bam"
        shell:
            "( [ -f {output} ] || ln -s -r {input} {output} ) "

# TODO: Make optional
rule conversionRate:
    input:
        "QC_metrics/{sample}.CHH.Mbias.txt"
    output:
        "QC_metrics/{sample}.conv.rate.txt"
    threads: 1
    shell: """
        awk '{{if(NR>1) {{M+=$4; UM+=$5}}}}END{{printf("{wildcards.sample}\\t%f\\n", 100*(1.0-M/(M+UM)))}}' {input} > {output}
        """


### bwameth ##########################################
if paired:
    rule bwameth:
        input:
            r1=fastq_dir + "/{sample}" + reads[0] + ".fastq.gz",
            r2=fastq_dir + "/{sample}" + reads[1] + ".fastq.gz"
        output:
            sbam=temp("bwameth/{sample}.sorted.bam")
        log:
            err="bwameth/logs/{sample}.map_reads.err",
            out="bwameth/logs/{sample}.map_reads.out"
        params:
            bwameth_index=bwameth_index
        threads: 20
        conda: CONDA_WGBS_ENV
        shell: """
            MYTEMP=$(mktemp -d "${{TMPDIR:-/tmp}}"/snakepipes.XXXXXXXXXX)
            bwameth.py --threads {threads} --reference "{params.bwameth_index}" "{input.r1}" "{input.r2}" 2> {log.err} | \
	        samtools sort -T "$MYTEMP"/{wildcards.sample} -m 3G -@ 4 -o "{output.sbam}"
            rm -rf "$MYTEMP"
            """
else:
    rule bwameth:
        input:
            r1=fastq_dir + "/{sample}" + reads[0] + ".fastq.gz",
        output:
            sbam=temp("bwameth/{sample}.sorted.bam")
        log:
            err="bwameth/logs/{sample}.map_reads.err",
            out="bwameth/logs/{sample}.map_reads.out"
        params:
            bwameth_index=bwameth_index
        threads: 20
        conda: CONDA_WGBS_ENV
        shell: """
            MYTEMP=$(mktemp -d "${{TMPDIR:-/tmp}}"/snakepipes.XXXXXXXXXX)
            bwameth.py --threads {threads} --reference "{params.bwameth_index}" "{input.r1}" 2> {log.err} | \
	        samtools sort -T "$MYTEMP/{wildcards.sample}" -m 3G -@ 4 -o "{output.sbam}"
            rm -rf "$MYTEMP"
            """


rule index_bam:
    input:
        "bwameth/{sample}.sorted.bam"
    output:
        temp("bwameth/{sample}.sorted.bam.bai")
    log:
        err="bwameth/logs/{sample}.index_bam.err",
        out="bwameth/logs/{sample}.index_bam.out"
    conda: CONDA_SHARED_ENV
    shell: """
        samtools index "{input}" >{log.out} 2>{log.err}
        """


rule markDupes:
    input:
        "bwameth/{sample}.sorted.bam",
        "bwameth/{sample}.sorted.bam.bai"
    output:
        "bwameth/{sample}.PCRrm.bam"
    log:
        err="bwameth/logs/{sample}.rm_dupes.err",
        out="bwameth/logs/{sample}.rm_dupes.out"
    threads: 10
    conda: CONDA_SAMBAMBA_ENV
    shell: """
        MYTEMP=$(mktemp -d "${{TMPDIR:-/tmp}}"/snakepipes.XXXXXXXXXX)
        sambamba markdup -t {threads} --tmpdir "$MYTEMP/{wildcards.sample}" "{input[0]}" "{output}" >{log.out} 2>{log.err}
        rm -rf "$MYTEMP"
        """


rule index_PCRrm_bam:
    input:
        "bwameth/{sample}.PCRrm.bam"
    output:
        "bwameth/{sample}.PCRrm.bam.bai"
    params:
    log:
        err="bwameth/logs/{sample}.index_PCRrm_bam.err",
        out="bwameth/logs/{sample}.index_PCRrm_bam.out"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: """
        samtools index "{input}" 1>{log.out} 2>{log.err}
        """


# TODO: I'm not sure how useful this really is. We could just run plotCoverage instead.
rule getRandomCpGs:
    output:
        temp("QC_metrics/randomCpG.bed")
    params:
        genome_fasta=genome_fasta
    run:
        import random
        random.seed(1234)  # Just to ensure reproducibility

        buf = []
        maxLen = 1000000  # In theory this could be changed
        chroms = []
        position = 0
        chars = 0
        lastChar = 'N'

        def addPosition(B, tid, pos, n):
            if len(B) < maxLen:
                B.append((tid, pos))
            else:
                x = random.randint(0, n)
                if x < maxLen:
                    B[x] = (tid, pos)

        for line in open(params['genome_fasta']):
            line = line.strip().split()[0]
            if line.startswith(">"):
                chroms.append(line[1:])
                lastChar = 'N'
                position = 0
                continue
            for c in line:
                if (lastChar == 'C' or lastChar == 'c') and (c == 'G' or c == 'g'):
                    addPosition(buf, len(chroms) - 1, position - 1, chars)
                    lastChar = 'N'
                    chars += 1
                else:
                    lastChar = c
                position += 1

        # Sort
        buf.sort()

        # write output
        if len(buf):
            o = open(output[0], "w")
            for tid, pos in buf:
                o.write("{}\t{}\t{}\n".format(chroms[tid], pos, pos + 2))
            o.close()


rule calc_Mbias:
    input:
        "bwameth/{sample}.PCRrm.bam",
        "bwameth/{sample}.PCRrm.bam.bai"
    output:
        "QC_metrics/{sample}.Mbias.txt"
    params:
        genome=genome_fasta
    log:
        out="QC_metrics/logs/{sample}.calc_Mbias.out"
    threads: 10
    conda: CONDA_WGBS_ENV
    shell: """
        MethylDackel mbias -@ {threads} {params.genome} {input[0]} QC_metrics/{wildcards.sample} 2> {output} > {log.out}
        """


rule calcCHHbias:
    input:
        "bwameth/{sample}.PCRrm.bam",
        "bwameth/{sample}.PCRrm.bam.bai"
    output:
        temp("QC_metrics/{sample}.CHH.Mbias.txt")
    params:
        genome=genome_fasta
    log:
        err="QC_metrics/logs/{sample}.calcCHHbias.err"
    threads: 10
    conda: CONDA_WGBS_ENV
    shell: """
        MethylDackel mbias -@ {threads} --CHH --noCpG --noSVG {params.genome} {input[0]} QC_metrics/{wildcards.sample} > {output} 2> {log.err}
        """


rule calc_GCbias:
    input:
        "bwameth/{sample}.PCRrm.bam",
        "bwameth/{sample}.PCRrm.bam.bai"
    output:
        GCbiasTXT="QC_metrics/{sample}.freq.txt",
        GCbiasPNG="QC_metrics/{sample}.GCbias." + plotFormat
    params:
        genomeSize=genome_size,
        twobitpath=genome_2bit,
        plotFormat=plotFormat
    log:
        out="QC_metrics/logs/{sample}.calc_GCbias.out"
    threads: 20
    conda: CONDA_SHARED_ENV
    shell: """
        computeGCBias -b {input[0]} --effectiveGenomeSize {params.genomeSize} -g {params.twobitpath} -l 300 --GCbiasFrequenciesFile {output.GCbiasTXT} -p {threads} --biasPlot {output.GCbiasPNG}
        """


# TODO: plotCoverage would be better
rule DepthOfCov:
    input:
        "bwameth/{sample}.PCRrm.bam",
        "bwameth/{sample}.PCRrm.bam.bai",
        "QC_metrics/randomCpG.bed"
    output:
        "QC_metrics/{sample}.doc.sample_summary",
    params:
        options="-F 'mapping_quality > 4 and not duplicate and not failed_quality_control' -c 0 -q 10",
        thresholds="-T 0 -T 1 -T 2 -T 5 -T 10 -T 15 -T 20 -T 30 -T 50"
    threads: 10
    log:
        err="QC_metrics/logs/{sample}.DepthOfCov.err",
    conda: CONDA_SAMBAMBA_ENV
    shell: """
        # sambamba depth returns a table of per-position depth metrics
        # The awk command keeps the header and averages these to global metrics
        sambamba depth region --combined {params.options} {params.thresholds} -t {threads} -m -L {input[2]} {input[0]} 2> {log.err} \
            | cut -f 4,6- \
            | awk '{{ if(NR==1) {{ \
                        print \
                      }} else {{ \
                        for(i=1; i<=NF; i++) vals[i]+=$i \
                      }} \
                   }} \
                   END{{ \
                     printf("%f", vals[1] / (NR - 1)); \
                     for(i=2; i<=NF; i++) {{ \
                       printf("\t%f", vals[i] / (NR - 1)) \
                     }} \
                     OFS=""; print("") \
                   }}' \
            > {output}
        """


rule DepthOfCovGenome:
    input:
        BAMS=expand("bwameth/{sample}.PCRrm.bam", sample=samples),
        BAIS=expand("bwameth/{sample}.PCRrm.bam.bai", sample=samples)
    output:
        "QC_metrics/genomeCoverage.sample_summary",
        "QC_metrics/genomeCoverage.png"
    threads: 20
    log:
        err="QC_metrics/logs/DepthOfCovGenome.err"
    conda: CONDA_SHARED_ENV
    shell: """
        plotCoverage -b {input.BAMS} -p {threads} --minMappingQuality 10 --smartLabels -o {output[1]} > {output[0]} 2> {log.err}
        """


rule get_flagstat:
    input:
        "bwameth/{sample}.PCRrm.bam"
    output:
        "QC_metrics/{sample}.flagstat"
    log:
        err="QC_metrics/logs/{sample}.get_flagstat.err"
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: "samtools flagstat {input} > {output} 2>{log.err}"


# TODO: CpG coverage stuff isn't very useful, use plotCoverage instead
rule produce_report:
    input:
        Coverage=calc_doc(skipDOC),
        ConversionRate=expand("QC_metrics/{sample}.conv.rate.txt", sample=samples),
        mbiasTXT=expand("QC_metrics/{sample}.Mbias.txt", sample=samples),
        fstat=expand("QC_metrics/{sample}.flagstat", sample=samples)
    output:
        QCrep='QC_metrics/QC_report.html'
    params:
        auxdir=os.path.join(outdir, "QC_metrics")
    conda: CONDA_RMD_ENV
    script: "../rscripts/WGBS_QC_report_template.Rmd"


if not noAutoMethylationBias:
    rule methyl_extract:
        input:
            "bwameth/{sample}.PCRrm.bam",
            "bwameth/{sample}.PCRrm.bam.bai",
            "QC_metrics/{sample}.Mbias.txt"
        output:
            "MethylDackel/{sample}_CpG.bedGraph"
        params:
            genome=genome_fasta,
            MethylDackelOptions=MethylDackelOptions
        log:
            err="MethylDackel/logs/{sample}.methyl_extract.err",
            out="MethylDackel/logs/{sample}.methyl_extract.out"
        threads: 10
        conda: CONDA_WGBS_ENV
        shell: """
            mi=$(cat {input[2]} | sed 's/Suggested inclusion options: //' )
            MethylDackel extract -o MethylDackel/{wildcards.sample} {params.MethylDackelOptions} $mi -@ {threads} {params.genome} {input[0]} 1>{log.out} 2>{log.err}
            """
else:
    rule methyl_extract:
        input:
            "bwameth/{sample}.PCRrm.bam",
            "bwameth/{sample}.PCRrm.bam.bai"
        output:
            "MethylDackel/{sample}_CpG.bedGraph"
        params:
            genome=genome_fasta,
            MethylDackelOptions=MethylDackelOptions
        log:
            err="MethylDackel/logs/{sample}.methyl_extract.err",
            out="MethylDackel/logs/{sample}.methyl_extract.out"
        threads: 10
        conda: CONDA_WGBS_ENV
        shell: """
            MethylDackel extract -o MethylDackel/{wildcards.sample} {params.MethylDackelOptions} -@ {threads} {params.genome} {input[0]} 1>{log.out} 2>{log.err}
            """


if blacklist is None:
    rule CpG_filt:
        input:
            "MethylDackel/{sample}_CpG.bedGraph"
        output:
            temp("MethylDackel/{sample}.CpG.filt2.bed")
        threads: 1
        conda: CONDA_WGBS_ENV
        shell: """
            echo -e "chr\tstart\tend\tBeta\tM\tU\tCov\tms" > {output}
            awk '{{if(NR>1) print $0, $5+$6, $1"_"$2}}' {input} | tr " " "\t" >> {output}
            """
else:
    rule CpG_filt:
        input:
            "MethylDackel/{sample}_CpG.bedGraph",
            blacklist
        output:
            temp("MethylDackel/{sample}.CpG.filt2.bed"),
            temp("MethylDackel/{sample}.CpG.filt2.bed.temp")
        log:
            err="MethylDackel/logs/{sample}.CpG_filt.err",
            out="MethylDackel/logs/{sample}.CpG_filt.out"
        threads: 1
        conda: CONDA_WGBS_ENV
        shell: """
            echo -e "chr\tstart\tend\tBeta\tM\tU\tCov\tms" > {output[1]}
            awk '{{if(NR>1) print $0, $5+$6, $1"_"$2}}' {input[0]} | tr " " "\t" >> {output[1]}
            bedtools intersect -v -a {output[1]} -b {input[1]} > {output[0]} 1>{log.out} 2>{log.err}
            """


# TODO: this is really slow, it shouldn't take more than a few minutes.
rule prep_for_stats:
    input:
        expressionFiles=expand("MethylDackel/{sample}.CpG.filt2.bed", sample=samples)
    output:
        MetileneIN='{}/metilene.IN.txt'.format(get_outdir("metilene"))
    params:
        sampleSheet=sampleSheet,
        importfunc="WGBSstats_functions.R"
    log:
        err='{}/logs/prep_for_stats.err'.format(get_outdir("metilene")),
    threads: 1
    conda: CONDA_WGBS_ENV
    script: "../rscripts/WGBS_mergeStats.R"


# TODO: allow changing smoothing parameters
# This is currently allotted 5GB per thread
rule DSS:
    input:
        bedGraphs=expand("MethylDackel/{sample}.CpG.filt2.bed", sample=samples)
    output:
        '{}/Stats_report.html'.format(get_outdir("DSS"))
    params:
        odir=get_outdir("DSS"),
        sampleSheet=sampleSheet,
        groups=metileneGroups,
        maxDist=maxDist,
        minCpGs=minCpGs,
        minMethDiff=minMethDiff,
        FDR=FDR
    threads: 10
    benchmark: '{}/.benchmark/DSS.benchmark'.format(get_outdir("DSS"))
    conda: CONDA_WGBS_ENV
    script: "../rscripts/WGBS_DSS.Rmd"


# metileneGroups is provided by the calling snakeFile
# These are NOT filtered
rule run_metilene:
    input:
        MetIN='{}/metilene.IN.txt'.format(get_outdir("metilene"))
    output:
        MetBed='{}/DMRs.txt'.format(get_outdir("metilene"))
    params:
        groups=metileneGroups,
        maxDist=maxDist,
        minCpGs=minCpGs,
        minMethDiff=minMethDiff,
        FDR=FDR
    log:
        err="{}/logs/run_metilene.err".format(get_outdir("metilene"))
    threads: 10
    benchmark: '{}/.benchmark/run_metilene.benchmark'.format(get_outdir("metilene"))
    conda: CONDA_WGBS_ENV
    shell: """
        echo -e "chrom\tstart\tend\tq-value\tmean methylation difference\tnCpGs\tp (MWU)\tp (2D KS)\tmean_{params.groups[0]}\tmean_{params.groups[1]}" > {output}
        metilene --groupA {params.groups[0]} \
                 --groupB {params.groups[1]} \
                 --maxdist {params.maxDist} \
                 --mincpgs {params.minCpGs} \
                 --minMethDiff {params.minMethDiff} \
                 --threads {threads} \
                 {input.MetIN} 2>{log.err} \
            | sort -k 1,1 -k2,2n >> {output.MetBed}
        """


# Annotates the DMRs and produces QC plots
rule metileneReport:
    input:
        '{}/DMRs.txt'.format(get_outdir("metilene")),
    output:
        HTML='{}/Stats_report.html'.format(get_outdir("metilene"))
    params:
        genes_gtf=genes_gtf,
        outdir=get_outdir("metilene"),
        sampleSheet=sampleSheet,
        minMethDiff=minMethDiff,
        FDR=FDR
    threads: 1
    benchmark: '{}/.benchmark/metileneReport.benchmark'.format(get_outdir("metilene"))
    conda: CONDA_WGBS_ENV
    script: "../rscripts/WGBS_metileneQC.Rmd"


rule bedGraphToBigWig:
    input: 
        "MethylDackel/{sample}_CpG.bedGraph",
        genome_index
    output:
        "MethylDackel/{sample}_CpG.methylation.bw",
        "MethylDackel/{sample}_CpG.coverage.bw"
    log:
        err='MethylDackel/logs/{sample}_bedGraphToBigWig.stderr'
    threads: 1
    conda: CONDA_SHARED_ENV
    shell: os.path.join(workflow_tools, "bedGraphToBigwig") + " {input[0]} {input[1]} {output[0]} {output[1]} 2> {log.err}"
