import subprocess as sp
import pytest
from ruamel.yaml import YAML
import pandas as pd

GENOME = "ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz"
GTF = "ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz"
RMSK = "http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz"
SPIKEINGENOME = "ftp://ftp.ensembl.org/pub/release-79/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa.gz"
SPIKEINGTF = "ftp://ftp.ensembl.org/pub/release-96/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.22.96.gtf.gz"
SMKOPTS = " --dryrun "

def parseSpOut(_s) -> int:
    '''
    parse subprocess run output.
    Take stdout, split, and take the 'jobnumber' field.
    The jobnumber field is assumed to be the list entry after 'total'
    Returns as int.
    '''
    for i in range(len(_s.stdout.split())-1):
        if 'total' in _s.stdout.split()[i]:
            return(int(_s.stdout.split()[i+1]))
    return (0)


def createTestData(fp, samples=9) -> None:
    yaml = YAML()
    yaml.boolean_representation = ['False', 'True']
    # single end folder
    (fp / 'SE').mkdir()
    # paired end folder
    (fp / 'PE').mkdir()
    # bam input folder
    (fp / 'bam_input' / 'filtered_bam').mkdir(parents=True)
    (fp / 'bam_input' / 'Sambamba').mkdir(parents=True)
    (fp / 'bam_input' / 'deepTools_qc' / 'bamPEFragmentSize').mkdir(parents=True)
    (fp / 'bam_input' / 'deepTools_qc' / 'bamPEFragmentSize' / 'fragmentSize.metric.tsv' ).touch()
    (fp / 'bam_input' / 'bamCoverage').mkdir(parents=True)
    # allelic bam input folder
    (fp / 'allelic_bam_input' / 'filtered_bam').mkdir(parents=True)
    (fp / 'allelic_bam_input' / 'allelic_bams').mkdir(parents=True)
    (fp / 'allelic_bam_input' / 'deepTools_qc' / 'bamPEFragmentSize').mkdir(parents=True)
    (fp / 'allelic_bam_input' / 'deepTools_qc' / 'bamPEFragmentSize' / 'fragmentSize.metric.tsv' ).touch()
    (fp / 'allelic_bam_input' / 'Sambamba').mkdir(parents=True)
    (fp / 'allelic_bam_input' / 'bamCoverage' / 'allele_specific').mkdir(parents=True)

    (fp / 'ref').mkdir()
    (fp / 'ref' / 'genes.gtf').touch()
    with open(fp / 'ref' / 'genome.fa', 'w') as f:
        f.write('>1\n')
        f.write('AAAAA\n')
        f.write('>2_spikein\n')
        f.write('TTTTT\n')
    with open(fp / 'ref' / 'genome.fa.fai', 'w') as f:
        f.write('1\5\n')
        f.write('2_spikein\t5\n')
    (fp / 'ref' / 'rmsk.txt').touch()
    (fp / 'ref' / 'genes.bed').touch()
    (fp / 'ref' / 'genes.slop.gtf').touch()
    (fp / 'ref' / 'spikein_genes.gtf').touch()
    (fp / 'ref' / 'genome.2bit').touch()
    (fp / 'ref' / 'splicesites.txt').touch()
    (fp / 'ref' / 'rar.bed').touch()
    (fp / 'ref' / 'decoys.txt').touch()
    (fp / 'ref' / 'cDNA_introns.joint.t2g').touch()
    (fp / 'ref' / 'genes.t2g').touch()

    (fp / 'allelic_input'/ 'Ngenome').mkdir(parents=True)
    (fp / 'allelic_input'/ 'file.vcf.gz').touch()
    (fp / 'allelic_input'/ 'snpfile.txt').touch()

    # samples
    for s in range(samples):
        sample = s+1
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
        (fp / "allelic_bam_input" / "allelic_bams" / "sample{}.allele_flagged.sorted.bam".format(sample)).touch()
        (fp / "allelic_bam_input" / "allelic_bams" / "sample{}.allele_flagged.sorted.bam.bai".format(sample)).touch()
        (fp / "allelic_bam_input" / "allelic_bams" / "sample{}.unassigned.sorted.bam".format(sample)).touch()
        (fp / "allelic_bam_input" / "allelic_bams" / "sample{}.unassigned.sorted.bam.bai".format(sample)).touch()
        (fp / "allelic_bam_input" / "filtered_bam" / "sample{}.filtered.bam".format(sample)).touch()
        (fp / "allelic_bam_input" / "filtered_bam" / "sample{}.filtered.bam.bai".format(sample)).touch()
        (fp / "allelic_bam_input" / "Sambamba" / "sample{}.markdup.txt".format(sample)).touch()
        (fp / "allelic_bam_input" / "bamCoverage" / "allele_specific" / "sample{}.genome1.seq_depth_norm.bw".format(sample)).touch()

    # Create organism.yaml
    orgyaml = {
        'genome_size': 2652783500,
        'genome_fasta': (fp / 'ref' / 'genome.fa').as_posix(),
        'genome_index': (fp / 'ref' / 'genome.fa.fai').as_posix(),
        'genome_2bit': (fp / 'ref' / 'genome.2bit').as_posix(),
        'bowtie2_index': (fp / 'ref' / 'genome').as_posix(),
        'hisat2_index': (fp / 'ref' / 'genome').as_posix(),
        'bwa_index': (fp / 'ref' / 'genome.fa').as_posix(),
        'bwa_mem2_index': (fp / 'ref' / 'genome.fa').as_posix(),
        'bwameth_index': (fp / 'ref' / 'genome.fa').as_posix(),
        'bwameth2_index': (fp / 'ref' / 'genome.fa').as_posix(),
        'known_splicesites': (fp / 'ref' / 'splicesites.txt').as_posix(),
        'star_index': (fp / 'ref').as_posix(),
        'salmon_index': (fp / 'ref').as_posix(),
        'salmon_velocity_index': (fp /  'ref').as_posix(),
        't2g_velocity': (fp / 'ref' / 'cDNA_introns.joint.t2g').as_posix(),
        'genes_bed': (fp / 'ref' / 'genes.bed').as_posix(),
        'genes_gtf': (fp / 'ref' / 'genes.gtf').as_posix(),
        'genes_t2g': (fp / 'ref' / 'genes.t2g').as_posix(),
        'spikein_genes_gtf' : (fp / 'ref' / 'spikein_genes.gtf').as_posix(),
        'extended_coding_regions_gtf': (fp / 'ref' / 'genes.slop.gtf').as_posix(),
        'blacklist_bed': (fp / 'ref' / 'rar.bed').as_posix(),
        'spikein_blacklist_bed': "",
        'ignoreForNormalization': "MT X Y",
        'rmsk_file': (fp / 'ref' / 'rmsk.txt').as_posix()
    }
    with open(fp / 'org.yaml', 'w') as of:
        yaml.dump(orgyaml, of)

    # create test samplesheet
    pd.DataFrame(
        [
            ['sample1', 'Control'],
            ['sample2', 'Control'],
            ['sample4', 'Treatment'],
            ['sample5', 'Treatment']
        ],
        columns = ['name', 'condition']
    ).to_csv(fp / 'sampleSheet.tsv', sep='\t', index=False)

    # create multicomp samplesheet
    pd.DataFrame(
        [
            ['sample1', 'Control', 'All'],
            ['sample2', 'Control', 'All'],
            ['sample3', 'Treatment', 'Group1'],
            ['sample4', 'Treatment', 'Group1'],
            ['sample5', 'Treatment', 'Group2'],
            ['sample6', 'Treatment', 'Group2'],
        ],
        columns = ['name', 'condition', 'group']
    ).to_csv(fp / 'sampleSheet_mc.tsv', sep='\t', index=False)
    # ChIP sample_config
    chip_dict ={
        'chip_dict': {
            'sample1': {'control': 'sample7', 'broad': False},
            'sample2': {'control': 'sample7', 'broad': False},
            'sample3': {'control': 'sample8', 'broad': False},
            'sample4': {'control': 'sample8', 'broad': False},
            'sample5': {'control': 'sample9', 'broad': False},
            'sample6': {'control': 'sample9', 'broad': False}
        }
    }
    with open(fp / 'chipdict.yaml', 'w') as f:
        yaml.dump(chip_dict, f)

    chip_dict_broad ={
        'chip_dict': {
            'sample1': {'control': 'sample7', 'broad': True},
            'sample2': {'control': 'sample7', 'broad': True},
            'sample3': {'control': 'sample8', 'broad': True},
            'sample4': {'control': 'sample8', 'broad': True},
            'sample5': {'control': 'sample9', 'broad': True},
            'sample6': {'control': 'sample9', 'broad': True}
        }
    }
    with open(fp / 'chipdict_broad.yaml', 'w') as f:
        yaml.dump(chip_dict_broad, f)

    chip_dict_no_control ={
        'chip_dict': {
            'sample1': {'control': None, 'broad': False},
            'sample2': {'control': None, 'broad': False},
            'sample3': {'control': None, 'broad': False},
            'sample4': {'control': None, 'broad': False},
            'sample5': {'control': None, 'broad': False},
            'sample6': {'control': None, 'broad': False}
        }
    }
    with open(fp / 'chipdict_noControl.yaml', 'w') as f:
        yaml.dump(chip_dict_no_control, f)


    chip_dict_simple ={
        'chip_dict': {
            'sample1': {'control': 'sample3', 'broad': False},
            'sample2': {'control': 'sample3', 'broad': False},
            'sample4': {'control': 'sample6', 'broad': False},
            'sample5': {'control': 'sample6', 'broad': False}
        }
    }
    with open(fp / 'chipdict_simple.yaml', 'w') as f:
        yaml.dump(chip_dict_simple, f)

@pytest.fixture(scope="session")
def ifs(tmp_path_factory):
    fp = tmp_path_factory.mktemp("data")
    createTestData(fp)
    return fp

class TestmakePairs():
    def test_default(self, ifs):
        ci = [
            "makePairs",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'output',
            ifs / 'org.yaml',
            '--VCFfile',
            ifs / 'allelic_input' / 'file.vcf.gz',
            '--strains',
            'strain1,strain2',
            '--snakemakeOptions',
            SMKOPTS
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 176

    def test_dag(self, ifs):
        ci = [
            "makePairs",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'output',
            '--DAG',
            ifs / 'org.yaml',
            '--VCFfile',
            ifs / 'allelic_input' / 'file.vcf.gz',
            '--strains',
            'strain1,strain2',
            '--snakemakeOptions',
            SMKOPTS
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 176

class TestCreateindices:
    def test_default(self, ifs):
        ci = [
            'createIndices',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            '--genome',
            GENOME,
            '--gtf',
            GTF,
            'genome'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 21
    def test_rmsk(self, ifs):
        ci = [
            'createIndices',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            '--genome',
            GENOME,
            '--gtf',
            GTF,
            'genome',
            '--rmskURL',
            RMSK,
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 22
    def test_DAG(self, ifs):
        ci = [
            'createIndices',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            '--genome',
            GENOME,
            '--gtf',
            GTF,
            'genome',
            '--rmskURL',
            RMSK,
            '--DAG'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 22
    def test_spikein(self, ifs):
        ci = [
            'createIndices',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            '--genome',
            GENOME,
            '--gtf',
            GTF,
            'genome',
            '--rmskURL',
            RMSK,
            '--DAG',
            '--spikeinGtfURL',
            SPIKEINGTF,
            '--spikeinGenomeURL',
            SPIKEINGENOME
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 27

class TestDNAmapping():
    def test_default(self, ifs):
        # PE
        ci = [
            "DNAmapping",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 143
    def test_properPairs(self, ifs):
        ci = [
            "DNAmapping",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--trim',
            '--mapq',
            '20',
            '--dedup',
            '--properPairs'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 152
    def test_bcExtract(self, ifs):
        ci = [
            "DNAmapping",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--trim',
            '--mapq',
            '20',
            '--dedup',
            '--properPairs',
            '--bcExtract'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 143
    def test_UMIDedup(self, ifs):
        ci = [
            "DNAmapping",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--trim',
            '--mapq',
            '20',
            '--UMIDedup',
            '--properPairs'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 161
    def test_UMIDedupbcExtract(self, ifs):
        ci = [
            "DNAmapping",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--trim',
            '--mapq',
            '20',
            '--UMIDedup',
            '--properPairs',
            '--bcExtract'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 152
    def test_DAG(self, ifs):
        ci = [
            "DNAmapping",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--DAG',
            '--trim',
            '--mapq',
            '20',
            '--UMIDedup',
            '--properPairs',
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 161
    def test_bwa(self, ifs):
        ci = [
            "DNAmapping",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--DAG',
            '--trim',
            '--aligner',
            'bwa'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 143
    def test_bwa2(self, ifs):
        ci = [
            "DNAmapping",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--DAG',
            '--trim',
            '--aligner',
            'bwa-mem2'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 143
    def test_se(self, ifs):
        ci = [
            "DNAmapping",
            '-i',
            ifs / 'SE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 125
    def test_seproperPairs(self, ifs):
        ci = [
            "DNAmapping",
            '-i',
            ifs / 'SE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--trim',
            '--mapq',
            '20',
            '--dedup',
            '--properPairs'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 134

class TestChIPseq:
    def test_default(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 59
    def test_nosamplesheet(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 37
    def test_nosamplesheet_genrich(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--snakemakeOptions',
            SMKOPTS,
            '--peakCaller',
            'Genrich',
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 35
    def test_broad(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict_broad.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 77
    def test_genrich(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            '--peakCaller',
            'Genrich',
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 56
    def test_seacr(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            '--peakCaller',
            'SEACR',
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 89
    def test_seacr_spikein(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            '--peakCaller',
            'SEACR',
            '--useSpikeInForNorm',
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 130
    def test_SE(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            '--singleEnd',
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 59
    def test_l2ratio(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            '--bigWigType',
            'log2ratio',
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 53
    def test_default_noInput(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict_noControl.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 38
    def test_genrich_noInput(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            '--peakCaller',
            'Genrich',
            ifs / 'org.yaml',
            ifs / 'chipdict_noControl.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 32
    def test_seacr_noInput(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            '--peakCaller',
            'SEACR',
            ifs / 'org.yaml',
            ifs / 'chipdict_noControl.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 62
    def test_seacr_spikein_noInput(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            '--peakCaller',
            'SEACR',
            '--useSpikeInForNorm',
            ifs / 'org.yaml',
            ifs / 'chipdict_noControl.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 83
    def test_frombam(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'outdir',
            '--fromBAM',
            ifs / 'bam_input' / 'filtered_bam',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 119
    def test_frombam_noInput(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'outdir',
            '--fromBAM',
            ifs / 'bam_input' / 'filtered_bam',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict_noControl.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 80
    def test_spikein(self, ifs):
        ci = [
            "ChIPseq",
            '--useSpikeInForNorm',
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 106
    def test_spikein_noInput(self, ifs):
        ci = [
            "ChIPseq",
            '--useSpikeInForNorm',
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict_noControl.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 65
    def test_spikeinfrombam(self, ifs):
        ci = [
            "ChIPseq",
            '--useSpikeInForNorm',
            '-d',
            ifs / 'outdir',
            '--fromBAM',
            ifs / 'bam_input' / 'filtered_bam',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 142
    def test_spikeinfrombamTSSnorm(self, ifs):
        ci = [
            "ChIPseq",
            '--useSpikeInForNorm',
            '--getSizeFactorsFrom',
            'TSS',
            '-d',
            ifs / 'outdir',
            '--fromBAM',
            ifs / 'bam_input' / 'filtered_bam',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 118
    def test_spikeinfrombaminputnorm(self, ifs):
        ci = [
            "ChIPseq",
            '--useSpikeInForNorm',
            '--getSizeFactorsFrom',
            'input',
            '-d',
            ifs / 'outdir',
            '--fromBAM',
            ifs / 'bam_input' / 'filtered_bam',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 115
    def test_allelic(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'allelic_bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict_simple.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 37
    def test_multicomp(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 81
    def test_multicomp_genrich(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            '--peakCaller',
            'Genrich',
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 82
    def test_multicomp_broad(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict_broad.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 99
    def test_multicomp_fromBam(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'outdir',
            '--fromBAM',
            ifs / 'bam_input' / 'filtered_bam',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 141
    def test_multicomp_fromBam_Genrich(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'outdir',
            '--fromBAM',
            ifs / 'bam_input' / 'filtered_bam',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            '--peakCaller',
            'Genrich',
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 142
    def test_multicomp_spikein(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--useSpikeInForNorm',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 126
    def test_multicomp_spikein_genrich(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--useSpikeInForNorm',
            '--snakemakeOptions',
            SMKOPTS,
            '--peakCaller',
            'Genrich',
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 127
    def test_multicomp_spikein_noInput(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--useSpikeInForNorm',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict_noControl.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 81
    def test_multicomp_spikein_noInput_Genrich(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--useSpikeInForNorm',
            '--snakemakeOptions',
            SMKOPTS,
            '--peakCaller',
            'Genrich',
            ifs / 'org.yaml',
            ifs / 'chipdict_noControl.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 79
    def test_multicomp_spikein_fromBam(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'outdir',
            '--fromBAM',
            ifs / 'bam_input' / 'filtered_bam',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--useSpikeInForNorm',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 162
    def test_multicomp_spikein_fromBam_genrich(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'outdir',
            '--fromBAM',
            ifs / 'bam_input' / 'filtered_bam',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--useSpikeInForNorm',
            '--snakemakeOptions',
            SMKOPTS,
            '--peakCaller',
            'Genrich',
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 163
    def test_multicomp_spikein_fromBam_noInput(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'outdir',
            '--fromBAM',
            ifs / 'bam_input' / 'filtered_bam',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--useSpikeInForNorm',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict_noControl.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 105
    def test_multicomp_spikein_fromBam_noInput_genrich(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'outdir',
            '--fromBAM',
            ifs / 'bam_input' / 'filtered_bam',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--useSpikeInForNorm',
            '--snakemakeOptions',
            SMKOPTS,
            '--peakCaller',
            'Genrich',
            ifs / 'org.yaml',
            ifs / 'chipdict_noControl.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 103
    def test_multicomp_fromBam_noInput_SEACR(self, ifs):
        ci = [
            "ChIPseq",
            '-d',
            ifs / 'outdir',
            '--fromBAM',
            ifs / 'bam_input' / 'filtered_bam',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            '--peakCaller',
            'SEACR',
            ifs / 'org.yaml',
            ifs / 'chipdict_noControl.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 120

class TestmRNAseq:
    def test_default(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 167
    def test_DE(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 168
    def test_rMats(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--rMats',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 170
    def test_almode(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '-m',
            'alignment'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 125
    def test_trim(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--trim'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 177
    def test_alfreemode(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '-m',
            'alignment-free,deepTools_qc'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 188
    def test_bcExtract(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--bcExtract',
            '--trim'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 168
    def test_bcExtractUMIdedup(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--bcExtract',
            '--trim',
            '--UMIDedup'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 177
    def test_multicomp(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '-m',
            'alignment,alignment-free',
            '--rMats',
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 165
    def test_SE(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'SE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 149
    def test_SEalmode(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'SE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '-m',
            'alignment'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 107
    def test_SEtrim(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'SE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--trim'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 158
    def test_SEalfreemode(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'SE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '-m',
            'alignment-free,deepTools_qc'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 169
    def test_SEfastqc(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'SE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--trim',
            '--fastqc'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 176
    def test_SEfrombam(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'bam_input' / 'filtered_bam',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--fromBAM'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 114
    def test_threeprime(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '-m',
            'three-prime-seq'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 189
    def test_threeprimeqc(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '-m',
            'three-prime-seq,deepTools_qc'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 232
    def test_allelic(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '-m',
            'allelic-mapping,deepTools_qc',
            '--VCFfile',
            ifs / 'allelic_input' / 'file.vcf.gz',
            '--strains',
            'strain1,strain2'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 253
    def test_allelicfrombam(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'allelic_bam_input' / 'filtered_bam',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--fromBAM',
            '-m',
            'allelic-mapping,deepTools_qc',
            '--SNPfile',
            ifs / 'allelic_input' / 'snpfile.txt',
            '--NMaskedIndex',
            ifs / 'allelic_input' / 'Ngenome'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 197
    def test_allelicDE(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '-m',
            'allelic-mapping,deepTools_qc',
            '--VCFfile',
            ifs / 'allelic_input' / 'file.vcf.gz',
            '--strains',
            'strain1,strain2'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 254
    def test_allelicDE_SNPfile(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '-m',
            'allelic-mapping,deepTools_qc',
            '--SNPfile',
            ifs / 'allelic_input' / 'snpfile.txt',
            '--NMaskedIndex',
            ifs / 'allelic_input' / 'Ngenome'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 252
    def test_allelicDEsinglestrain(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '-m',
            'allelic-mapping,deepTools_qc',
            '--VCFfile',
            ifs / 'allelic_input' / 'file.vcf.gz',
            '--strains',
            'strain1'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 254
    def test_allelicDEalfree(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--VCFfile',
            ifs / 'allelic_input' / 'file.vcf.gz',
            '--strains',
            'strain1',
            '-m',
            'allelic-mapping,deepTools_qc,alignment-free'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 330
    def test_allelic_count_fromBam_singlecomp(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'allelic_bam_input' / 'allelic_bams',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            '--fromBAM',
            '--bamExt',
            '.sorted.bam',
            ifs / 'org.yaml',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '-m',
            'allelic-counting'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 105
    def test_allelic_count_fromBam_multicomp(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'allelic_bam_input' / 'allelic_bams',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            '--fromBAM',
            '--bamExt',
            '.sorted.bam',
            ifs / 'org.yaml',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '-m',
            'allelic-counting'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 108
    def test_allelic_mapping_fromBam_multicomp(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'allelic_bam_input' / 'filtered_bam',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            '--fromBAM',
            '--bamExt',
            '.filtered.bam',
            '-m',
            'allelic-mapping,deepTools_qc',
            ifs / 'org.yaml',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--SNPfile',
            ifs / 'allelic_input' / 'snpfile.txt',
            '--NMaskedIndex',
            ifs / 'allelic_input' / 'Ngenome'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 201
    def test_allelic_alfree_multicomp(self, ifs):
        ci = [
            "mRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '-m',
            'allelic-mapping,deepTools_qc,alignment-free',
            '--VCFfile',
            ifs / 'allelic_input' / 'file.vcf.gz',
            '--strains',
            'strain1'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 331

class TestncRNAseq():
    def test_default(self, ifs):
        ci = [
            "ncRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 137
    def test_DE(self, ifs):
        ci = [
            "ncRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 139
    def test_SE(self, ifs):
        ci = [
            "ncRNAseq",
            '-i',
            ifs / 'SE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 120
    def test_frombam(self, ifs):
        ci = [
            "ncRNAseq",
            '-i',
            ifs / 'bam_input' / 'filtered_bam',
            '--fromBAM',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 94
    def test_multicomp(self, ifs):
        ci = [
            "ncRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 142

class TestscRNAseq():
    def test_default(self, ifs):
        ci = [
            "scRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--mode',
            'STARsolo',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 178
    def test_skipvelo(self, ifs):
        ci = [
            "scRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--mode',
            'STARsolo',
            '--skipVelocyto',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 159
    def test_alevin(self, ifs):
        ci = [
            "scRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--mode',
            'Alevin',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 78
    def test_alevinskipvelo(self, ifs):
        ci = [
            "scRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--mode',
            'Alevin',
            '--skipVelocyto',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 66

class TestWGBS():
    def test_default(self, ifs):
        ci = [
            "WGBS",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 154
    def test_no_sampleSheet(self, ifs):
        ci = [
            "WGBS",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 150
    def test_bwameth2(self, ifs):
        ci = [
            "WGBS",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--aligner',
            'bwameth2',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 154
    def test_trimgcbias(self, ifs):
        ci = [
            "WGBS",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'outdir',
            '--trim',
            '--GCbias',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 155
    def test_frombam(self, ifs):
        ci = [
            "WGBS",
            '-i',
            ifs / 'bam_input' / 'filtered_bam',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--fromBAM',
            '--GCbias',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 110
    def test_frombamfqc(self, ifs):
        ci = [
            "WGBS",
            '-i',
            ifs / 'bam_input' / 'filtered_bam',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--fromBAM',
            '--GCbias',
            '--fastqc',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 110
    def test_frombamskipqc(self, ifs):
        ci = [
            "WGBS",
            '-i',
            ifs / 'bam_input' / 'filtered_bam',
            '-o',
            ifs / 'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--fromBAM',
            '--skipBamQC',
            '--fastqc',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 59

class TestATAC():
    def test_default(self, ifs):
        ci = [
            "ATACseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 63
    def test_no_sampleSheet(self, ifs):
        ci = [
            "ATACseq",
            '-d',
            ifs / 'bam_input',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 47
    def test_genrich(self, ifs):
        ci = [
            "ATACseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--peakCaller',
            'Genrich',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 74
    def test_HMMRATAC(self, ifs):
        ci = [
            "ATACseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--peakCaller',
            'HMMRATAC',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 73
    def test_sieve(self, ifs):
        ci = [
            "ATACseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--maxFragmentSize',
            '120',
            '--qval',
            '0.1',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 63
    def test_frombam(self, ifs):
        ci = [
            "ATACseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            '--fromBAM',
            ifs / 'bam_input' / 'filtered_bam',
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 114
    def test_multicomp_default(self, ifs):
        ci = [
            "ATACseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 79
    def test_multicomp_genrich(self, ifs):
        ci = [
            "ATACseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--peakCaller',
            'Genrich',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 92
    def test_multicomp_HMMRATAC(self, ifs):
        ci = [
            "ATACseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--peakCaller',
            'HMMRATAC',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 89
    def test_multicomp_sieve(self, ifs):
        ci = [
            "ATACseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--maxFragmentSize',
            '120',
            '--qval',
            '0.1',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 79
    def test_multicomp_frombam(self, ifs):
        ci = [
            "ATACseq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            '--fromBAM',
            ifs / 'bam_input' / 'filtered_bam',
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 130

class TestHIC():
    def test_default(self, ifs):
        ci = [
            "HiC",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'output',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 102

    def test_ice(self, ifs):
        ci = [
            "HiC",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'output',
            '--snakemakeOptions',
            SMKOPTS,
            '--correctionMethod',
            'ICE',
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 111

    def test_trim(self, ifs):
        ci = [
            "HiC",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'output',
            '--snakemakeOptions',
            SMKOPTS,
            '--trim',
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 111

    def test_dpnii(self, ifs):
        ci = [
            "HiC",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'output',
            '--snakemakeOptions',
            SMKOPTS,
            '--enzyme',
            'DpnII',
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 102

    def test_notad(self, ifs):
        ci = [
            "HiC",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'output',
            '--snakemakeOptions',
            SMKOPTS,
            '--noTAD',
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 93

    def test_bwamem2(self, ifs):
        ci = [
            "HiC",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'output',
            '--snakemakeOptions',
            SMKOPTS,
            '--aligner',
            'bwa-mem2',
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 102

class Testpreprocessing():
    def test_default(self, ifs):
        ci = [
            "preprocessing",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'output',
            '--snakemakeOptions',
            SMKOPTS,
            '--fastqc',
            '--optDedupDist',
            '2500'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 84

    def test_DAG(self, ifs):
        ci = [
            "preprocessing",
            '-i',
            ifs / 'PE',
            '-o',
            ifs / 'output',
            '--snakemakeOptions',
            SMKOPTS,
            '--fastqc',
            '--optDedupDist',
            '2500',
            '--DAG'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 84
