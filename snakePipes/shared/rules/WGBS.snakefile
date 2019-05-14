import os
import re
from operator import is_not
import tempfile

# N.B., this needs a `get_outdir()` function from the importing Snakefile

###symlink bams if this is the starting point
if fromBam:
    rule link_bam:
        input:
            indir + "/{sample}" + bam_ext
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
            bwameth.py --threads {threads} --reference {params.bwameth_index} {input.r1} {input.r2} 2> {log.err} | \
	        samtools sort -T ${{TMPDIR}}/{wildcards.sample} -m 3G -@ 4 -o {output.sbam}
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
            bwameth.py --threads {threads} --reference {params.bwameth_index} {input.r1} 2> {log.err} | \
	        samtools sort -T ${{TMPDIR}}/{wildcards.sample} -m 3G -@ 4 -o {output.sbam}
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
    shell: "samtools index {input} 1>{log.out} 2>{log.err}"


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
    shell: "sambamba markdup -t {threads} --tmpdir ${{TMPDIR}}{wildcards.sample} {input[0]} {output} 1>{log.out} 2>{log.err}"


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
    shell: "samtools index {input} 1>{log.out} 2>{log.err}"


# TODO: I'm not sure how useful this really is. We could just run plotCoverage instead.
rule getRandomCpGs:
    output:
        temp("aux_files/randomCpG.bed")
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
        GCbiasPNG="QC_metrics/{sample}.GCbias." + plot_format
    params:
        genomeSize=genome_size,
        twobitpath=genome_2bit,
        plot_format=plot_format
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
        intList if intList else "aux_files/randomCpG.bed"
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
        # The awk command keeps the header and averages these to global metric
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
        auxdir=os.path.join(outdir, "aux_files")
    conda: CONDA_RMD_ENV
    script: "../rscripts/WGBS_QC_report_template.Rmd"


if mbias_ignore=="auto":
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
            MethylDackelOptions=MethylDackelOptions,
            mbias_ignore=mbias_ignore
        log:
            err="MethylDackel/logs/{sample}.methyl_extract.err",
            out="MethylDackel/logs/{sample}.methyl_extract.out"
        threads: 10
        conda: CONDA_WGBS_ENV
        shell: """
            MethylDackel extract -o MethylDackel/{wildcards.sample} {params.mbias_ignore} -@ {threads} {params.genome} {input[0]} 1>{log.out} 2>{log.err}
            """


if blackList is None:
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
            methTab="MethylDackel/{sample}_CpG.bedGraph",
            blackListF=blackList
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


rule prep_for_stats:
    input:
        expressionFiles=expand("MethylDackel/{sample}.CpG.filt2.bed", sample=samples)
    output:
        LimmaData=temp('{}/limdat.LG.RData'.format(get_outdir("limma"))),
        MetileneIN=temp('{}/metilene.IN.txt'.format(get_outdir("metilene")))
    params:
        sampleSheet=sampleSheet,
        importfunc="WGBSstats_functions.R"
    log:
        err='{}/logs/prep_for_stats.err'.format(get_outdir("limma")),
    threads: 1
    conda: CONDA_WGBS_ENV
    script: "../rscripts/WGBS_mergeStats.R"


# TODO, make temp?
# TODO: integrate
#rule CpG_stats:
#    input: 
#        Limdat='{}/limdat.LG.RData'.format(get_outdir("limma"))
#    output:
#        RDatAll='{}/singleCpG.RData'.format(get_outdir("limma")),
#        sinfo='{}/sessionInfo.txt'.format(get_outdir("limma"))
#    params:
#        statdir=os.path.join(outdir,'{}'.format(get_outdir("limma"))),
#        datdir=os.path.join(outdir,'{}'.format(get_outdir("limma"))),
#        script=os.path.join(workflow_rscripts,'WGBSpipe.singleCpGstats.limma.R'),
#        sampleSheet=sampleSheet,
#        diff=minAbsDiff,
#        fdr=FDR,
#        importfunc = os.path.join(workflow_rscripts, "WGBSstats_functions.R")
#    log:
#        err='{}/logs/CpG_stats.err'.format(get_outdir("limma")),
#        out='{}/logs/CpG_stats.out'.format(get_outdir("limma"))
#    threads: 1
#    conda: CONDA_WGBS_ENV
#    shell: """
#        Rscript --no-save --no-restore {params.script} {params.statdir} {params.sampleSheet} {params.datdir} {params.diff} {params.fdr} {params.importfunc} 1>{log.out} 2>{log.err}
#        """


# TODO: make temp?
# TODO: integrate
rule CpG_report:
    input: 
        Limdat='{}/limdat.LG.RData'.format(get_outdir("limma")),
        sinfo='{}/singleCpG.RData'.format(get_outdir("limma"))
    output:
        html='{}/Stats_report.html'.format(get_outdir("limma"))
    params:
        statdir=os.path.join(outdir,'{}'.format(get_outdir("limma"))),
        sampleSheet=sampleSheet,
        importfunc = os.path.join(workflow_rscripts, "WGBSstats_functions.R"),
        stat_cat="single_CpGs",
        rmd_in=os.path.join(workflow_rscripts,"WGBS_stats_report_template.Rmd"),
        rmd_out=os.path.join(outdir,"aux_files", "WGBS_stats_report_template.Rmd"),
        outFull=lambda wildcards,output: os.path.join(outdir,output.html)
    log:
        err='{}/logs/stats_report.err'.format(get_outdir("limma")),
        out='{}/logs/stats_report.out'.format(get_outdir("limma"))
    conda: CONDA_RMD_ENV
    threads: 1
    shell: """
        cp -v {params.rmd_in} {params.rmd_out}
        Rscript -e 'rmarkdown::render("{params.rmd_out}", params=list(outdir="{params.statdir}", input_func="{params.importfunc}", stat_category="{params.stat_cat}",sample_sheet="{params.sampleSheet}"), output_file="{params.outFull}")' 1>{log.out} 2>{log.err}
        """


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
    conda: CONDA_WGBS_ENV
    script: "../rscripts/WGBS_metileneQC.Rmd"


# TODO: imdF needs to be relabeled
# TODO: rewrite shell script
#rule get_CG_per_int:
#    input:
#        intList=intList,
#        refG=genome_fasta,
#        imdF="aux_files/"+re.sub('.fa*','.CpG.bed',os.path.basename(genome_fasta))
#    output:
#        outList=run_int_aggStats(intList,False)
#    log:
#        err="aux_files/logs/get_CG_per_int.err"
#    params:
#        auxshell=lambda wildcards,input,output: ';'.join(["bedtools intersect -wa -a "+ input.imdF + " -b " + bli + ' > ' + oli  for bli,oli in zip(input.intList,output.outList) ])
#    threads: 1
#    conda: CONDA_WGBS_ENV
#    shell: "{params.auxshell} 2>{log.err}"


## TODO: rewrite
#rule intAgg_stats:
#    input:
#        Limdat='{}/limdat.LG.RData'.format(get_outdir("limma")),
#        intList=intList,
#        refG=genome_fasta,
#        sampleSheet=sampleSheet,
#        aux=run_int_aggStats(intList,False)
#    output:
#        outFiles=run_int_aggStats(intList,sampleSheet),
#        sinfo='{}/sessionInfo.txt'.format(get_outdir("aggregate_stats_limma"))
#    params:
#        auxshell=lambda wildcards,input:';'.join(['Rscript --no-save --no-restore ' + os.path.join(workflow_rscripts,'WGBSpipe.interval_stats.limma.R ') + os.path.join(outdir,'{}'.format(get_outdir("aggregate_stats_limma"))) + ' ' + li +' '+ aui +' ' + os.path.join(outdir,input.Limdat) + ' '  + input.sampleSheet + ' ' + str(minAbsDiff) + ' ' + str(FDR) + ' ' + os.path.join(workflow_rscripts, "WGBSstats_functions.R")  for li,aui in zip(intList,[os.path.join(outdir,"aux_files",re.sub('.fa',re.sub('.bed','.CpGlist.bed',os.path.basename(x)),os.path.basename(genome_fasta))) for x in intList])])
#    log:
#        err="{}/logs/intAgg_stats.err".format(get_outdir("aggregate_stats_limma")),
#        out="{}/logs/intAgg_stats.out".format(get_outdir("aggregate_stats_limma"))
#    threads: 1
#    conda: CONDA_WGBS_ENV
#    shell: "{params.auxshell} 1>{log.out} 2>{log.err}"


## TODO:
#rule intAgg_report:
#    input: 
#        outFiles='{}/sessionInfo.txt'.format(get_outdir("aggregate_stats_limma"))
#    output:
#        html='{}/Stats_report.html'.format(get_outdir("aggregate_stats_limma"))
#    params:
#        statdir=os.path.join(outdir,'{}'.format(get_outdir("aggregate_stats_limma"))),
#        sampleSheet=sampleSheet,
#        importfunc = os.path.join(workflow_rscripts, "WGBSstats_functions.R"),
#        stat_cat="user_intervals",
#        rmd_in=os.path.join(workflow_rscripts,"WGBS_stats_report_template.Rmd"),
#        rmd_out=os.path.join(outdir,"aux_files", "WGBS_stats_report_template.Rmd"),
#        outFull=lambda wildcards,output: os.path.join(outdir,output.html)
#    log:
#        err='{}/logs/stats_report.err'.format(get_outdir("aggregate_stats_limma")),
#        out='{}/logs/stats_report.out'.format(get_outdir("aggregate_stats_limma"))
#    conda: CONDA_RMD_ENV
#    threads: 1
#    shell: "cp -v {params.rmd_in} {params.rmd_out} ;Rscript -e 'rmarkdown::render(\"{params.rmd_out}\", params=list(outdir=\"{params.statdir}\", input_func=\"{params.importfunc}\", stat_category=\"{params.stat_cat}\",sample_sheet=\"{params.sampleSheet}\"), output_file=\"{params.outFull}\")' 1>{log.out} 2>{log.err}"


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
