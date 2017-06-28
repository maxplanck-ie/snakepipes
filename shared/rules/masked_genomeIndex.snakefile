
from os.path import join, dirname


GENOMEDIR = os.path.dirname(genome_fasta)
BASENAME = genome
print(GENOMEDIR)
print(BASENAME)
# define snpgenome_dir
if allele_hybrid == "dual":
    SNPdir = "snp_genome/" + strains[0] + "_" + \
                    strains[1] + "_dual_hybrid.based_on_" + \
                    BASENAME + "_N-masked"
else:
    SNPdir = "snp_genome/" + strains[0] + "_" + "_N-masked"

if allele_hybrid == 'dual':
    rule snpgenome:
        input:
            genome = GENOMEDIR,
            snps = VCFfile
        output:
            genome1 = "snp_genome/" + strains[0] + '_SNP_filtering_report.txt',
            genome2 = "snp_genome/" + strains[1] + '_SNP_filtering_report.txt',
            snpgenome_dir = "snp_genome/" + strains[0] + "_" + \
                            strains[1] + "_dual_hybrid.based_on_" + \
                            BASENAME + "_N-masked"
        params:
            strain1 = strains[0],
            strain2 = strains[1]
        log: "snp_genome/SNPsplit_createSNPgenome.log"
        shell:
            SNPsplit_path +"SNPsplit_genome_preparation"
            " --dual_hybrid --genome_build {BASENAME} "
            " --reference_genome {input.genome} --vcf_file {input.snps}"
            " --strain {params.strain1} --strain2 {params.strain2}"
else:
    rule snpgenome:
        input:
            genome = GENOMEDIR,
            snps = VCFfile
        output:
            genome1 = "snp_genome/" + strains[0] + '_SNP_filtering_report.txt',
            snpgenome_dir = "snp_genome/" + strains[0] + "_" + "_N-masked"
        params:
            strain1 = strains[0]
        log: "snp_genome/SNPsplit_createSNPgenome.log"
        shell:
            SNPsplit_path +"SNPsplit_genome_preparation"
            " --genome_build {BASENAME} "
            " --reference_genome {input.genome} --vcf_file {input.snps}"
            " --strain {params.strain1}"

if mapping_prg == "STAR":
    rule star_index:
        input:
            genome = expand(SNPdir + '/{file}', file = INFILES)
        output:
            index = "snp_genome/star_Nmasked/Genome"
        log:
            "snp_genome/star_Nmasked/star.index.log"
        threads:
            10
        params:
            gtf=genes_gtf
        run:
            # Write stderr and stdout to the log file.
            shell("/package/STAR-2.5.2b/bin/STAR"
                  " --runThreadN {threads}"
                  " --runMode genomeGenerate"
                  " --genomeDir " + "star_Nmasked"
                  " --genomeFastaFiles {input.genome}"
                  " --sjdbGTFfile {params.gtf}"
                  " > {log} 2>&1")
else:
    print("Only STAR is implemented for allele-specific mapping")
