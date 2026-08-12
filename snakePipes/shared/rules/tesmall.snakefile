import os
import glob
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
import seaborn as sns
import textwrap

# TEsmall writes ~20 files per sample straight into TEsmallOut/ with no
# subdirectories -- this sorts them by type after a run so the directory is
# actually readable. Final deliverables (count_summary*, report.html,
# TEsmall.done) are left at the top level; everything else is per-sample
# intermediate/diagnostic output.
TESMALL_OUTPUT_CATEGORIES = {
    "log": (".log",),
    "bam": (".bam",),
    "anno": (".anno", ".anno.rlen.info"),
    "rinfo": (".rinfo",),
    "bedgraph": (".bedgraph",),
    "cca_fa": ("_cca.fa", ".unaligned.cca.fa"),
    "fastq": (".fastq",),
}


def organize_tesmall_output(teout_dir):
    for category in TESMALL_OUTPUT_CATEGORIES:
        os.makedirs(os.path.join(teout_dir, category), exist_ok=True)

    for fname in os.listdir(teout_dir):
        fpath = os.path.join(teout_dir, fname)
        if not os.path.isfile(fpath):
            continue
        for category, suffixes in TESMALL_OUTPUT_CATEGORIES.items():
            if fname.endswith(suffixes):
                shutil.move(fpath, os.path.join(teout_dir, category, fname))
                break

Adapterseq = "../../workflows/smRNAseq/reads_adapters_set.fasta"

# tesmall_db (--dbfolder) and tesmall_genome (-g) both come straight from the
# organism YAML createIndices wrote (see workflows/createIndices/Snakefile) --
# no need to guess the build from the --genome label here.

def get_samples(Samples):
    all_fastq = []
    for file in Samples:
        FqPaths = os.path.join(outdir, file)
        print(f"Generated file path: {FqPaths}")  # Debugging line
        all_fastq.append(FqPaths)
    return all_fastq

def justify_text(text, line_width):
    words = text.split( )
    lines = textwrap.wrap(text, width=line_width)
    justified_lines = []
    
    for line in lines:
        words_in_line = line.split()
        if len(words_in_line) > 1:
            spaces_needed = line_width - sum(len(w) for w in words_in_line)
            space_between = spaces_needed // (len(words_in_line) - 1)
            extra_spaces = spaces_needed % (len(words_in_line) - 1)
            
            justified_line = ""
            for i, word in enumerate(words_in_line):
                justified_line += word
                if i < len(words_in_line) - 1:
                    justified_line += " " * (space_between + (1 if i < extra_spaces else 0))
            justified_lines.append(justified_line)
        else:
            justified_lines.append(line)
    
    return "\n".join(justified_lines)

def Plotting(OutFile):
    countFile = os.path.join(outdir, 'TEsmallOut', 'count_summary.txt')
    df = pd.read_csv(countFile, sep='\t')
    samples = df.columns[2:]
    

    result_list = []

    for sample in samples:
        sample_data = {
            'Sample': sample,
            'anti_TE': df[df['ftype'] == 'anti_TE'][sample].sum(),
            'exon': df[df['ftype'] == 'exon'][sample].sum(),
            'hairpin': df[df['ftype'] == 'hairpin'][sample].sum(),
            'intron': df[df['ftype'] == 'intron'][sample].sum(),
            'miRNA': df[df['ftype'] == 'miRNA'][sample].sum(), 
            'piRNA_cluster': df[df['ftype'] == 'piRNA_cluster'][sample].sum(),
            'sense_TE': df[df['ftype'] == 'sense_TE'][sample].sum(),
            'structural_RNA': df[df['ftype'] == 'structural_RNA'][sample].sum(),
        }
        result_list.append(sample_data)

    result = pd.DataFrame(result_list)

    result.set_index('Sample', inplace=True)

    total = result.sum(axis=1)

    fraction_df = result.div(total, axis=0)

    plt.figure(figsize=(12, 8))
    ax = fraction_df.plot(kind='bar', stacked=True, figsize=(12, 8))

    for index, value in enumerate(total):
        ax.text(index, 1.02, f'{value/1e6:.1f}M', ha='center', va='bottom')
    
    x_labels = ax.get_xticklabels()
    wrapped_labels = [textwrap.fill(label.get_text(), width=15) for label in x_labels]
    
    ax.set_xticklabels(wrapped_labels)

    ax.set_ylabel('Fraction', fontsize=14)
    ax.set_xlabel('Sample', fontsize=14)
    ax.set_ylim(0, 1.2)

    ax.legend(title='Categories', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
    plt.xticks(rotation=0, ha='center', fontsize=12)
    plt.tight_layout()
    plt.savefig(OutFile)

if pairedEnd:
    rule fastp:
        input:
            fastq_indir_trim+"/{sample}"+reads[0]+".fastq.gz",
            fastq_indir_trim+"/{sample}"+reads[1]+".fastq.gz"
        output:
            #"FASTQ_fastp/{sample}"+".fastq.gz",
            os.path.join(outdir,"FASTQ_fastp","{sample}.fastq.gz"),
            "FASTQ_fastp/{sample}fastp.json",
            "FASTQ_fastp/{sample}fastp.html"
        benchmark:
            "FASTQ_fastp/.benchmark/fastp.{sample}.benchmark"
        threads: lambda wildcards: 8 if 8<max_thread else max_thread
        conda: CONDA_SHARED_ENV
        shell:"""
            fastp -w {threads} -i "{input[0]}" -I "{input[1]}"  -m --merged_out "{output[0]}" --include_unmerged -j "{output[1]}" -h "{output[2]}"
        """
else:
    rule fastp:
        input:
            fastq_indir_trim+"/{sample}"+reads[0]+".fastq.gz"
        output:
            os.path.join(outdir,"FASTQ_fastp","{sample}.fastq.gz"),
            #"FASTQ_fastp/{sample}"+".fastq.gz",
            "FASTQ_fastp/{sample}fastp.json",
            "FASTQ_fastp/{sample}fastp.html"
            
        benchmark:
            "FASTQ_fastp/.benchmark/fastp.{sample}.benchmark"
        threads: lambda wildcards: 8 if 8<max_thread else max_thread
        conda: CONDA_SHARED_ENV
        shell: """
            fastp -w {threads} -i "{input[0]}" --adapter_fasta {Adapterseq}  -o "{output[0]}" -j "{output[1]}" -h "{output[2]}"
            """

rule Tesmall_run:
    input:
        fqIn = expand(os.path.join(outdir,"FASTQ_fastp","{sample}.fastq.gz"), sample=samples)
    output:
        os.path.join(outdir,'TEsmallOut','TEsmall.done')
    threads: 16
    params:
        db = tesmall_db,
        genome_version = tesmall_genome,
        minlen = tesmallMinLen,
        maxlen = tesmallMaxLen,
        extra_opts = str(tesmallOptions or ''),
        outputdir = os.path.join(outdir, 'TEsmallOut')
    conda: CONDA_SMRNA_ENV
    shell:'''
        cd {params.outputdir}
        echo "TEsmall -f {input.fqIn} --dbfolder {params.db} -g {params.genome_version} -m {params.minlen} -M {params.maxlen} -p {threads} {params.extra_opts}"
        TEsmall -f {input.fqIn} --dbfolder {params.db} -g {params.genome_version} -m {params.minlen} -M {params.maxlen} -p {threads} {params.extra_opts}
        touch {output}
    '''

rule Plot_smRNA_stats:
    input:
        os.path.join(outdir, 'TEsmallOut', 'TEsmall.done')
    output:
        plot = os.path.join(outdir, 'TEsmallOut', 'count_summary_plot.pdf'),
        # organize_tesmall_output (below) sorts TEsmall's raw per-sample output into
        # TEsmallOut/bam/ by moving files there -- declared here so downstream rules
        # (deepTools QC) can depend on the genome-aligned bam without a checkpoint.
        genome_bam = expand(os.path.join(outdir, 'TEsmallOut', 'bam', '{sample}.genome.bam'), sample=samples)
    run:
        Plotting(output.plot)
        organize_tesmall_output(os.path.join(outdir, 'TEsmallOut'))

# TEsmall's own genome bam isn't sorted -- deepTools needs it sorted, and the shared
# deepTools_RNA.snakefile rules all read from filtered_bam/{sample}.filtered.bam, the
# convention every other workflow's aligner output is normalized to. Indexing is
# handled by the generic samtools_index_filtered rule (umi_tools.snakefile).
rule tesmall_genome_bam_filtered:
    input:
        bam = os.path.join(outdir, 'TEsmallOut', 'bam', '{sample}.genome.bam')
    output:
        bam = "filtered_bam/{sample}.filtered.bam"
    benchmark:
        "filtered_bam/.benchmark/tesmall_genome_bam_filtered.{sample}.benchmark"
    threads: lambda wildcards: 4 if 4<max_thread else max_thread
    conda: CONDA_SHARED_ENV
    shell: "samtools sort -@ {threads} -o {output.bam} {input.bam}"