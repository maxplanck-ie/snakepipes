import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
import seaborn as sns
import textwrap

Adapterseq = "/data/manke/processing1/navandar/Softwares/snakepipes/snakePipes/workflows/smRNAseq/reads_adapters_set.fasta"
Organisms = {

    'mm10_gencodeM19':'mm10',
    'dm6':'dm6',
    'mm39_ens106':'mm39',
    'hg38':'hg38',
    'GRCz11':'GRCz11'
}

org = next((value for key, value in Organisms.items() if key == genome), genome)

def get_organism(organism):
    DB_selection = tesmall_db
    return DB_selection

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
    # Plotting the count summary:
    countFile = os.path.join(outdir, 'TEsmallOut', 'count_summary.txt')
    df = pd.read_csv(countFile, sep='\t')
    samples = df.columns[2:]
    
    # Initialize an empty list to collect the sample data
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

    # Convert the list of sample data into a DataFrame
    result = pd.DataFrame(result_list)

    # Set the Sample column as index
    result.set_index('Sample', inplace=True)

    # Calculate total for each sample
    total = result.sum(axis=1)

    # Calculate fractions
    fraction_df = result.div(total, axis=0)

    # Plotting
    plt.figure(figsize=(12, 8))
    ax = fraction_df.plot(kind='bar', stacked=True, figsize=(12, 8))

    # Add totals on top of each bar
    for index, value in enumerate(total):
        ax.text(index, 1.02, f'{value/1e6:.1f}M', ha='center', va='bottom')
    
    x_labels = ax.get_xticklabels()
    wrapped_labels = [textwrap.fill(label.get_text(), width=15) for label in x_labels]
    
    ax.set_xticklabels(wrapped_labels)

    ax.set_ylabel('Fraction', fontsize=14)
    ax.set_xlabel('Sample', fontsize=14)
    ax.set_ylim(0, 1.2)  # Adjust the y-limit to fit the totals on top

    ax.legend(title='Categories', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=12)
    plt.xticks(rotation=0, ha='center', fontsize=12)
    plt.tight_layout()

    # Save the plot
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
        #fqIn = expand("FASTQ_fastp/{sample}"+".fastq.gz", sample=samples)
        fqIn = expand(os.path.join(outdir,"FASTQ_fastp","{sample}.fastq.gz"), sample=samples)
    output:
        os.path.join(outdir,'TEsmallOut','TEsmall.done')
    threads: 16
    params:
        db = get_organism(organism = org),
        outputdir = os.path.join(outdir, 'TEsmallOut')
    conda: CONDA_SMRNA_ENV
    shell:'''
        cd {params.outputdir}
        echo "TEsmall -f {input.fqIn} --dbfolder {params.db} -g {org} -p {threads}"
        TEsmall -f {input.fqIn} --dbfolder {params.db} -g {org} -p {threads}
        touch {output}
    '''

rule Plot_smRNA_stats:
    input:
        os.path.join(outdir, 'TEsmallOut', 'TEsmall.done')
    output:
        os.path.join(outdir,'TEsmallOut','count_summary_plot.pdf')
    run:
        Plotting(output[0])