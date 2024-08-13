from os.path import join, dirname
import glob

GENOMEDIR = os.path.dirname(genome_fasta)
genome_alias = os.path.splitext(os.path.basename(genome))[0]
BASENAME = genome_alias
# define snpgenome_dir
if allele_hybrid == 'dual':
    SNPdir = "snp_genome/" + strains[0] + "_" + \
                    strains[1] + "_dual_hybrid.based_on_" + \
                    BASENAME + "_N-masked"
else:
    SNPdir = "snp_genome/" + strains[0] + "_N-masked"

def getref_fileList(dir):
    fl = glob.glob(dir + "/*.fa")
    flist = ','.join(fl)
    return(flist)


## Create masked genome
if allele_hybrid == 'dual':
    rule create_snpgenome:
        input:
            genome = GENOMEDIR
        output:
            genome1 = "snp_genome/" + strains[0] + '_SNP_filtering_report.txt',
            genome2 = "snp_genome/" + strains[1] + '_SNP_filtering_report.txt',
            snpgenome_dir = directory(SNPdir),
            SNPFile = SNPFile
        params:
            strain1 = strains[0],
            strain2 = strains[1],
            SNPpath = os.path.abspath(VCFfile)
        log:
            out = "SNPsplit_createSNPgenome.out",
            err = "SNPsplit_createSNPgenome.err"
        conda: CONDA_SHARED_ENV
        shell:
            " ( [ -d snp_genome ] || mkdir -p snp_genome ) && cd snp_genome &&"
            " SNPsplit_genome_preparation"
            " --dual_hybrid --genome_build {BASENAME}"
            " --reference_genome {input.genome} --vcf_file {params.SNPpath}"
            " --strain {params.strain1} --strain2 {params.strain2} > {log.out} 2> {log.err}"
            "&& cd ../"
else:
    rule create_snpgenome:
        input:
            genome = GENOMEDIR
        output:
            genome1 = "snp_genome/" + strains[0] + '_SNP_filtering_report.txt',
            snpgenome_dir = directory(SNPdir),
            SNPFile = SNPFile
        params:
            strain1 = strains[0],
            SNPpath = os.path.abspath(VCFfile),

            temp_out=temp("all_SNPs_" + strains[0] + "_" + genome_alias + ".txt.gz"),
            out_bname=os.path.basename(SNPFile)
        log:
            out = "SNPsplit_createSNPgenome.out",
            err = "SNPsplit_createSNPgenome.err"
        conda: CONDA_SHARED_ENV
        shell:
            " ( [ -d snp_genome ] || mkdir -p snp_genome ) && cd snp_genome &&"
            " SNPsplit_genome_preparation"
            " --genome_build {BASENAME}"
            " --reference_genome {input.genome} --vcf_file {params.SNPpath}"
            " --strain {params.strain1} > {log.out} 2> {log.err} "
            #&& cp "
            #"{params.temp_out} {params.out_bname} >> {log.out} 2>> {log.err} "
            "&& cd ../"

if aligner == "STAR":
    rule star_index:
        input:
            snpgenome_dir = SNPdir
        output:
            star_index_allelic
        log:
            out = "snp_genome/star_Nmasked/star.index.out",
            err = "snp_genome/star_Nmasked/star.index.err"
        threads:
            10
        params:
            gtf=genes_gtf
        conda: CONDA_RNASEQ_ENV
        shell:
            "STAR"
            " --runThreadN {threads}"
            " --runMode genomeGenerate"
            " --genomeDir " + "snp_genome/star_Nmasked"
            " --genomeFastaFiles {input.snpgenome_dir}/*.fa"
            " --sjdbGTFfile {params.gtf}"
            " > {log.out} 2> {log.err}"

elif aligner == "Bowtie2":
    rule bowtie2_index:
        input:
            snpgenome_dir = SNPdir
        output:
            bowtie2_index_allelic
        log:
            out = "snp_genome/bowtie2_Nmasked/bowtie2.index.out",
            err = "snp_genome/bowtie2_Nmasked/bowtie2.index.err"
        threads: lambda wildcards: 10 if 10<max_thread else max_thread
        params:
            filelist = getref_fileList(SNPdir),
            idxbase = "snp_genome/bowtie2_Nmasked/Genome"
        conda: CONDA_DNA_MAPPING_ENV
        shell:
            "bowtie2-build"
            " --threads {threads}"
            " {params.filelist}"
            " {params.idxbase}"
            " > {log.out} 2> {log.err}"
else:
    print("Only STAR and Bowtie2 are implemented for allele-specific mapping")
