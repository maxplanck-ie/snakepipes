# bed_annot is only used in this file
# bed_filtered is used by deepTools

rule filter_gtf_spikein:
    input:
        gtf = spikein_genes_gtf
    output:
        gtf = "Annotation_spikein/genes.filtered.gtf"
    params:
        pattern = "" if not filterGTF else filterGTF
    shell: """
        if [ -z {params.pattern} ] ; then
            ln -s {input.gtf} {output.gtf}
        else
            grep {params.pattern} {input.gtf} > {output.gtf}
        fi
        """


# Given a GTF file, optionally filter it and produce the following files:
# 
# Annotation/genes.filtered.bed
#	A BED version of the filtered GTF with transcript entries. This is used by
#	deepTools
#
# Annotation/genes.filtered.t2g
#	Mapping of transcript to gene IDs, used by R
#
# Annotation/genes.filtered.symbol
#	Gene ID -> gene name mapping, used by R

rule spikein_gtf_to_files:
    input:
        gtf = "Annotation_spikein/genes.filtered.gtf"
    output:
        "Annotation_spikein/genes.filtered.t2g",
        "Annotation_spikein/genes.filtered.symbol",
        "Annotation_spikein/genes.filtered.bed"
    run:
        import shlex
        import re

        t2g = open(output[0], "w")
        symbol = open(output[1], "w")
        GTFdict = dict()

        for line in open(input.gtf):
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            annos = re.split(''';(?=(?:[^'"]|'[^']*'|"[^"]*")*$)''', cols[8]) 
            if cols[2] == "gene":
                # get the gene_name and gene_id values
                gene_id = None
                gene_name = None
                for anno in annos:
                    anno = shlex.split(anno.strip(), " ")
                    if len(anno) == 0:
                        continue
                    if anno[0] == "gene_id":
                        gene_id = anno[1]
                    elif anno[0] == "gene_name":
                        gene_name = anno[1]
                if gene_id:
                    symbol.write("{}\t{}\n".format(gene_id, "" if not gene_name else gene_name))
            elif cols[2] == "transcript" or 'RNA' in cols[2]:
                # get the gene_id and transcript_id values
                gene_id = None
                transcript_id = None
                gene_name = ""
                for anno in annos:
                    anno = shlex.split(anno.strip(), " ")
                    if len(anno) == 0:
                        continue
                    if anno[0] == "gene_id":
                        gene_id = anno[1]
                    elif anno[0] == "transcript_id":
                        transcript_id = anno[1]
                    elif anno[0] == "gene_name":
                        gene_name = anno[1]
                if transcript_id:
                    t2g.write("{}\t{}\t{}\n".format(transcript_id, "" if not gene_id else gene_id, gene_name))
                    # chrom, start, end, strand, exon width and exon start offset
                    GTFdict[transcript_id] = [cols[0], cols[3], cols[4], cols[6], [], []]
            elif cols[2] == "exon":
                # get the transcript_id
                transcript_id = None
                for anno in annos:
                    anno = shlex.split(anno.strip(), " ")
                    if len(anno) == 0:
                        continue
                    if anno[0] == "transcript_id":
                        transcript_id = anno[1]
                if transcript_id and transcript_id in GTFdict:
                    exonWidth = int(cols[4]) - int(cols[3]) + 1
                    exonOffset = int(cols[3]) - int(GTFdict[transcript_id][1])
                    GTFdict[transcript_id][4].append(str(exonWidth))
                    GTFdict[transcript_id][5].append(str(exonOffset))

        t2g.close()
        symbol.close()

        BED = open(output[2], "w")
        for k, v in GTFdict.items():
            # sort the starts and sizes together
            v[5] = [int(x) for x in v[5]]
            v[4] = [int(x) for x in v[4]]
            blockSizes = [str(x) for _,x in sorted(zip(v[5], v[4]))]
            blockStarts = sorted(v[5])
            blockStarts = [str(x) for x in blockStarts]
            BED.write("{}\t{}\t{}\t{}\t.\t{}\t{}\t{}\t255,0,0\t{}\t{}\t{}\n".format(v[0],  # chrom
                                                                               v[1],  # start
                                                                               v[2],  # end
                                                                               k,
                                                                               v[3],  # strand
                                                                               v[1],  # start
                                                                               v[2],  # end
                                                                               len(v[4]),  # blockCount
                                                                               ",".join(blockSizes),  # blockSizes
                                                                               ",".join(blockStarts)))  # blockStarts
        BED.close()


## make standard annotation


rule TSS_to_windows:
    input:
        bed = "Annotation_spikein/genes.filtered.bed",
        bam = expand("split_bam/{sample}_spikein.bam",sample=samples)[0]
    output:
        TSS_bed = "Annotation_spikein/TSS.filtered.bed"
    params:
        outfile = outdir+"/Annotation_spikein/TSS.filtered.bed",
        script = maindir+"/shared/tools/TSS_to_windows.py",
        size = 5000
    conda: CONDA_pysam_ENV
    shell: """
            python {params.script} -inf {input.bed} -outf {output.TSS_bed} -size {params.size} -bam {input.bam}
              """
        
