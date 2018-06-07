.. _WGBS:

WGBS
============

Input requirements and outputs:
-------------------------------------------
This pipeline requires paired-end reads fastq files and a bisulfite converted genome as inputs. Methylation values per CpG dinucleotide are extracted after read mapping. 
Statistical analysis of differential methylation is performed if user provides a bed file with genomic intervals to aggregate methylation values over. De novo DMR calling, statistical evaluation and annotation with nearest gene are performed automatically. Various quality metrics are collected and a QC report is output.
