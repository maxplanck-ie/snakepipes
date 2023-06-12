import subprocess as sp
import pytest

def parseSpOut(_s) -> int:
    '''
    parse subprocess run output.
    Take stdout, split, and take the 'jobnumber' field.
    Returns as int.
    snakemake's last line in a quiet dryrun has:
    totaljobs, jobcount, min threads and max threads.
    '''
    return (int(_s.stdout.split()[-3]))

def createTestData(fp, samples=6) -> None:
    (fp / 'SE').mkdir()
    (fp / 'PE').mkdir()
    (fp / 'bam_input' / 'filtered_bam').mkdir(parents=True)
    (fp / 'bam_input' / 'Sambamba').mkdir(parents=True)
    (fp / 'bam_input' / 'deeptools_qc' / 'bamPEFragmentSize').mkdir(parents=True)
    (fp / 'bam_input' / 'deeptools_qc' / 'bamPEFragmentSize' / 'fragmentSize.metric.tsv' ).touch()
    (fp / 'bam_input' / 'bamCoverage').mkdir(parents=True)
    (fp / 'allelic_bam_input' / 'filtered_bam').mkdir(parents=True)
    (fp / 'allelic_bam_input' / 'allelic_bams').mkdir(parents=True)
    (fp / 'allelic_bam_input' / 'deeptools_qc' / 'bamPEFragmentSize').mkdir(parents=True)
    (fp / 'allelic_bam_input' / 'deeptools_qc' / 'bamPEFragmentSize' / 'fragmentSize.metric.tsv' ).touch()
    (fp / 'allelic_bam_input' / 'Sambamba').mkdir(parents=True)
    (fp / 'allelic_bam_input' / 'bamCoverage' / 'allele_specific').mkdir(parents=True)
    
    (fp / 'ref').mkdir()
    (fp / 'ref' / 'genes.gtf').touch()
    (fp / 'ref' / 'genome.fa').touch()
    (fp / 'ref' / 'genome.fa.fai').touch()
    (fp / 'ref' / 'rmsk.txt').touch()
    (fp / 'ref' / 'genes.bed').touch()
    (fp / 'ref' / 'spikein_genes.gtf').touch()
    (fp / 'ref' / 'genome.2bit').touch()

    (fp / 'allelic_input'/ 'Ngenome').mkdir(parents=True)
    (fp / 'allelic_input'/ 'file.vcf.gz').touch()
    (fp / 'allelic_input'/ 'snpfile.txt').touch()

    # samples
    for sample in range(samples):
        # SE
        (fp / "SE" / "sample{}_R1.fastq.gz".format(sample)).touch()
        
        # PE
        (fp / "PE" / "sample{}_R1.fastq.gz".format(sample)).touch()
        (fp / "PE" / "sample{}_R2.fastq.gz".format(sample)).touch()
        
        # bam_input
        (fp / "bam_input" / "sample{}.bam".format(sample)).touch()
        (fp / "bam_input" / "filtered_bam" / "sample{}.filtered.bam".format(sample)).touch()
        (fp / "bam_input" / "filtered_bam" / "sample{}.filtered.bam.bai".format(sample)).touch()
        (fp / "bam_input" / "Sambamba" / "sample{}.markdup.txt".format(sample)).touch()
        (fp / "bam_input" / "bamCoverage" / "sample{}.filtered.seq_depth_norm.bw".format(sample)).touch()
        
        # allelic_bams
        (fp / "allelic_bam_input" / "sample{}.bam".format(sample)).touch()
        (fp / "allelic_bam_input" / "allelic_bams" / "sample{}.genome1.sorted.bam".format(sample)).touch()
        (fp / "allelic_bam_input" / "allelic_bams" / "sample{}.genome1.sorted.bam.bai".format(sample)).touch()
        (fp / "allelic_bam_input" / "allelic_bams" / "sample{}.genome2.sorted.bam".format(sample)).touch()
        (fp / "allelic_bam_input" / "allelic_bams" / "sample{}.genome2.sorted.bam.bai".format(sample)).touch()
        (fp / "allelic_bam_input" / "filtered_bam" / "sample{}.filtered.bam".format(sample)).touch()
        (fp / "allelic_bam_input" / "filtered_bam" / "sample{}.filtered.bam.bai".format(sample)).touch()
        (fp / "allelic_bam_input" / "Sambamba" / "sample{}.markdup.txt".format(sample)).touch()
        (fp / "allelic_bam_input" / "bamCoverage" / "allele_specific" / "sample{}.genome1.seq_depth_norm.bw".format(sample)).touch()



@pytest.fixture(scope="session")
def inputfiles(tmp_path_factory):
    fp = tmp_path_factory.mktemp("data")
    createTestData(fp)
    return fp

class TestCreateindices:

    GENOME = "ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz"
    GTF = "ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz"
    RMSK = "http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz"
    SPIKEINGENOME = "ftp://ftp.ensembl.org/pub/release-79/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa.gz"
    SPIKEINGTF = "ftp://ftp.ensembl.org/pub/release-96/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.22.96.gtf.gz"
    SMKOPTS = " --dryrun --conda-prefix /tmp -q "
    OUTPUT = "output"
    ORG = ".ci_stuff/organism.yaml"

    def test_createIndices(self):
        cic1 = [
            'createIndices',
            '-o',
            self.OUTPUT,
            '--snakemakeOptions',
            self.SMKOPTS,
            '--genome',
            self.GENOME,
            '--gtf',
            self.GTF,
            'genome'
        ]
        cic2 = cic1 + [
            '--rmskURL',
            self.RMSK,
        ]
        cic3 = cic2 + [
            '--DAG'
        ]
        cic4 = cic3 + [
            '--spikeinGtfURL',
            self.SPIKEINGTF,
            '--spikeinGenomeURL',
            self.SPIKEINGENOME
        ]

        s1 = sp.run(cic1, capture_output=True, text=True)
        s2 = sp.run(cic2, capture_output=True, text=True)
        s3 = sp.run(cic3, capture_output=True, text=True)
        s4 = sp.run(cic4, capture_output=True, text=True)

        assert (parseSpOut(s1), s1.returncode) == (17, 0)
        assert parseSpOut(s2) == 18
        assert s2.returncode == 0
        assert parseSpOut(s3) == 18
        assert s3.returncode == 0
        assert parseSpOut(s4) == 23
        assert s4.returncode == 0

    def test_dnamapping(self, inputfiles):
        # PE
        dnape1 = [
            "DNA-mapping",
            '-i',
            inputfiles / 'PE',
            '-o',
            self.OUTPUT,
            '--snakemakeOptions',
            self.SMKOPTS,
            self.ORG
        ]
        dnape2 = dnape1 + "--trim --mapq 20 --dedup --properPairs".split(' ')
        dnape3 = dnape1 +"--trim --mapq 20 --dedup --properPairs --bcExtract".split(' ')
        dnape4 = dnape1 + "--trim --mapq 20 --UMIDedup --properPairs --bcExtract".split(' ')
        dnape5 = dnape1 + "--trim --mapq 20 --UMIDedup --properPairs".split(' ')
        dnape6 = dnape1 + "--DAG --trim --mapq 20 --UMIDedup --properPairs".split(' ')
        dnape7 = dnape1 + "--DAG --trim --aligner bwa".split(' ')
        dnape8 = dnape1 + "--DAG --trim --aligner bwa-mem2".split(' ')

        s1 = sp.run(dnape1, capture_output=True, text=True)
        s2 = sp.run(dnape2, capture_output=True, text=True)
        s3 = sp.run(dnape3, capture_output=True, text=True)
        s4 = sp.run(dnape4, capture_output=True, text=True)
        s5 = sp.run(dnape5, capture_output=True, text=True)
        s6 = sp.run(dnape6, capture_output=True, text=True)
        s7 = sp.run(dnape7, capture_output=True, text=True)
        s8 = sp.run(dnape8, capture_output=True, text=True)
        
        # SE
        dnase1 = [
            'DNA-mapping',
            '-i',
            inputfiles / 'SE',
            '-o',
            self.OUTPUT,
            '--snakemakeOptions',
            self.SMKOPTS,
            self.ORG
        ]
        dnase2 = dnase1 + "--trim --mapq 20 --dedup --properPairs".split(' ')

        s9 = sp.run(dnase1, capture_output=True, text=True)
        s10 = sp.run(dnase2, capture_output=True, text=True)

        # PE
        assert (parseSpOut(s1), s1.returncode) == (98, 0)
        assert (parseSpOut(s2), s2.returncode) == (104, 0)
        assert (parseSpOut(s3), s3.returncode) == (98, 0)
        assert (parseSpOut(s4), s4.returncode) == (104, 0)
        assert (parseSpOut(s5), s5.returncode) == (110, 0)
        assert (parseSpOut(s6), s6.returncode) == (110, 0)
        assert (parseSpOut(s7), s7.returncode) == (98, 0)
        assert (parseSpOut(s8), s8.returncode) == (98, 0)
        # SE
        assert (parseSpOut(s9), s7.returncode) == (86, 0)
        assert (parseSpOut(s10), s8.returncode) == (92, 0)
