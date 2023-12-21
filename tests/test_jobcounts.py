import subprocess as sp
import pytest
import yaml
import pandas as pd


GENOME = "ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz"
GTF = "ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz"
RMSK = "http://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz"
SPIKEINGENOME = "ftp://ftp.ensembl.org/pub/release-79/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa.gz"
SPIKEINGTF = "ftp://ftp.ensembl.org/pub/release-96/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.22.96.gtf.gz"
SMKOPTS = " --dryrun --conda-prefix /tmp -q "


def parseSpOut(_s) -> int:
    '''
    parse subprocess run output.
    Take stdout, split, and take the 'jobnumber' field.
    Returns as int.
    snakemake's last line in a quiet dryrun has:
    totaljobs, jobcount, min threads and max threads.
    '''
    try:
        return (int(_s.stdout.split()[-3]))
    except IndexError:
        return (0)

def createTestData(fp, samples=9) -> None:
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
        yaml.dump(orgyaml, of, default_flow_style=False)

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
            'sample1': {'control': 'sample7', 'broad': 'False'},
            'sample2': {'control': 'sample7', 'broad': 'False'},
            'sample3': {'control': 'sample8', 'broad': 'False'},
            'sample4': {'control': 'sample8', 'broad': 'False'},
            'sample5': {'control': 'sample9', 'broad': 'False'},
            'sample6': {'control': 'sample9', 'broad': 'False'}
        }
    }
    with open(fp / 'chipdict.yaml', 'w') as f:
        yaml.dump(chip_dict, f, default_flow_style=False, default_style=None)

    chip_dict ={
        'chip_dict': {
            'sample1': {'control':, 'broad': 'False'},
            'sample2': {'control':, 'broad': 'False'},
            'sample3': {'control':, 'broad': 'False'},
            'sample4': {'control':, 'broad': 'False'},
            'sample5': {'control':, 'broad': 'False'},
            'sample6': {'control':, 'broad': 'False'}
        }
    }
    with open(fp / 'chipdict_noControl.yaml', 'w') as f:
        yaml.dump(chip_dict, f, default_flow_style=False, default_style=None)


    chip_dict ={
        'chip_dict': {
            'sample1': {'control': 'sample3', 'broad': 'False'},
            'sample2': {'control': 'sample3', 'broad': 'False'},
            'sample4': {'control': 'sample6', 'broad': 'False'},
            'sample5': {'control': 'sample6', 'broad': 'False'}
        }
    }
    with open(fp / 'chipdict_simple.yaml', 'w') as f:
        yaml.dump(chip_dict, f, default_flow_style=False, default_style=None)

@pytest.fixture(scope="session")
def ifs(tmp_path_factory):
    fp = tmp_path_factory.mktemp("data")
    createTestData(fp)
    return fp

class TestCreateindices:
    def test_default(self):
        ci = [
            'createIndices',
            '-o',
            'outdir',
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
    def test_rmsk(self):
        ci = [
            'createIndices',
            '-o',
            'outdir',
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
    def test_DAG(self):
        ci = [
            'createIndices',
            '-o',
            'outdir',
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
    def test_spikein(self):
        ci = [
            'createIndices',
            '-o',
            'outdir',
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
            "DNA-mapping",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
            "DNA-mapping",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
            "DNA-mapping",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
            "DNA-mapping",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
            "DNA-mapping",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
            "DNA-mapping",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
            "DNA-mapping",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
            "DNA-mapping",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
            "DNA-mapping",
            '-i',
            ifs / 'SE',
            '-o',
            'outdir',
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
            "DNA-mapping",
            '-i',
            ifs / 'SE',
            '-o',
            'outdir',
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
            "ChIP-seq",
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
        assert parseSpOut(_p) == 77
    def test_genrich(self, ifs):
        ci = [
            "ChIP-seq",
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
        assert parseSpOut(_p) == 74
    def test_SE(self, ifs):
        ci = [
            "ChIP-seq",
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
        assert parseSpOut(_p) == 77
    def test_l2ratio(self, ifs):
        ci = [
            "ChIP-seq",
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
        assert parseSpOut(_p) == 71
    def test_default_noInput(self, ifs):
        ci = [
            "ChIP-seq",
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
        assert parseSpOut(_p) == 60
    def test_genrich_noInput(self, ifs):
        ci = [
            "ChIP-seq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--peakCaller Genrich'
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict_noControl.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 60
    def test_frombam(self, ifs):
        ci = [
            "ChIP-seq",
            '-d',
            'outdir',
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
        assert parseSpOut(_p) == 102
    def test_frombam_noInput(self, ifs):
        ci = [
            "ChIP-seq",
            '-d',
            'outdir',
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
        assert parseSpOut(_p) == 102
    def test_spikein(self, ifs):
        ci = [
            "ChIP-seq",
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
        assert parseSpOut(_p) == 89
    def test_spikein_noInput(self, ifs):
        ci = [
            "ChIP-seq",
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
        assert parseSpOut(_p) == 89
    def test_spikeinfrombam(self, ifs):
        ci = [
            "ChIP-seq",
            '--useSpikeInForNorm',
            '-d',
            'outdir',
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
        assert parseSpOut(_p) == 113
    def test_spikeinfrombamTSSnorm(self, ifs):
        ci = [
            "ChIP-seq",
            '--useSpikeInForNorm',
            '--getSizeFactorsFrom',
            'TSS',
            '-d',
            'outdir',
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
        assert parseSpOut(_p) == 92
    def test_spikeinfrombaminputnorm(self, ifs):
        ci = [
            "ChIP-seq",
            '--useSpikeInForNorm',
            '--getSizeFactorsFrom',
            'input',
            '-d',
            'outdir',
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
        assert parseSpOut(_p) == 90
    def test_allelic(self, ifs):
        ci = [
            "ChIP-seq",
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
            "ChIP-seq",
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
        assert parseSpOut(_p) == 89
    def test_multicomp_genrich(self, ifs):
        ci = [
            "ChIP-seq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--peakCaller Genrich',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 89
    def test_multicomp_fromBam(self, ifs):
        ci = [
            "ChIP-seq",
            '-d',
            'outdir',
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
        assert parseSpOut(_p) == 89
    def test_multicomp_fromBam_Genrich(self, ifs):
        ci = [
            "ChIP-seq",
            '-d',
            'outdir',
            '--fromBAM',
            ifs / 'bam_input' / 'filtered_bam',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--peakCaller Genrich',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 89
    def test_multicomp_spikein(self, ifs):
        ci = [
            "ChIP-seq",
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
        assert parseSpOut(_p) == 89
    def test_multicomp_spikein_genrich(self, ifs):
        ci = [
            "ChIP-seq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--useSpikeInForNorm',
            '--peakCaller Genrich',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 89
    def test_multicomp_spikein_noInput(self, ifs):
        ci = [
            "ChIP-seq",
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
        assert parseSpOut(_p) == 89
    def test_multicomp_spikein_noInput_Genrich(self, ifs):
        ci = [
            "ChIP-seq",
            '-d',
            ifs / 'bam_input',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--useSpikeInForNorm',
            '--peakCaller Genrich',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict_noControl.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 89
    def test_multicomp_spikein_fromBam(self, ifs):
        ci = [
            "ChIP-seq",
            '-d',
            'outdir',
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
        assert parseSpOut(_p) == 89
    def test_multicomp_spikein_fromBam_genrich(self, ifs):
        ci = [
            "ChIP-seq",
            '-d',
            'outdir',
            '--fromBAM',
            ifs / 'bam_input' / 'filtered_bam',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--peakCaller Genrich',
            '--useSpikeInForNorm',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 89
    def test_multicomp_spikein_fromBam_noInput(self, ifs):
        ci = [
            "ChIP-seq",
            '-d',
            'outdir',
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
        assert parseSpOut(_p) == 89
    def test_multicomp_spikein_fromBam_noInput_genrich(self, ifs):
        ci = [
            "ChIP-seq",
            '-d',
            'outdir',
            '--fromBAM',
            ifs / 'bam_input' / 'filtered_bam',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--peakCaller Genrich',
            '--useSpikeInForNorm',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
            ifs / 'chipdict_noControl.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 89

class TestmRNAseq:
    def test_default(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 116
    def test_DE(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 117
    def test_rMats(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 119
    def test_almode(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 86
    def test_trim(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 123
    def test_alfreemode(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 131
    def test_bcExtract(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 117
    def test_bcExtractUMIdedup(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 123
    def test_multicomp(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 117
    def test_SE(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'SE',
            '-o',
            'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml',
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 104
    def test_SEalmode(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'SE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 74
    def test_SEtrim(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'SE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 110
    def test_SEalfreemode(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'SE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 118
    def test_SEfastqc(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'SE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 122
    def test_SEfrombam(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'bam_input' / 'filtered_bam',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 81
    def test_threeprime(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 132
    def test_threeprimeqc(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 163
    def test_allelic(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 175
    def test_allelicfrombam(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'allelic_bam_input' / 'filtered_bam',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 137
    def test_allelicDE(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 176
    def test_allelicDE_SNPfile(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 174
    def test_allelicDEsinglestrain(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 176
    def test_allelicDEalfree(self, ifs):
        ci = [
            "mRNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 228

class TestncRNAseq():
    def test_default(self, ifs):
        ci = [
            "noncoding-RNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 95
    def test_DE(self, ifs):
        ci = [
            "noncoding-RNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 97
    def test_SE(self, ifs):
        ci = [
            "noncoding-RNA-seq",
            '-i',
            ifs / 'SE',
            '-o',
            'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 84
    def test_frombam(self, ifs):
        ci = [
            "noncoding-RNA-seq",
            '-i',
            ifs / 'bam_input' / 'filtered_bam',
            '--fromBAM',
            '-o',
            'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 67
    def test_multicomp(self, ifs):
        ci = [
            "noncoding-RNA-seq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet_mc.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 100

class TestscRNAseq():
    def test_default(self, ifs):
        ci = [
            "scRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
            '--mode',
            'STARsolo',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 124
    def test_skipvelo(self, ifs):
        ci = [
            "scRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 111
    def test_alevin(self, ifs):
        ci = [
            "scRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
            '--mode',
            'Alevin',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 54
    def test_alevinskipvelo(self, ifs):
        ci = [
            "scRNAseq",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 45

class TestWGBS():
    def test_default(self, ifs):
        ci = [
            "WGBS",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
            '--sampleSheet',
            ifs / 'sampleSheet.tsv',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 106
    def test_bwameth2(self, ifs):
        ci = [
            "WGBS",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 106
    def test_trimgcbias(self, ifs):
        ci = [
            "WGBS",
            '-i',
            ifs / 'PE',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 107
    def test_frombam(self, ifs):
        ci = [
            "WGBS",
            '-i',
            ifs / 'bam_input' / 'filtered_bam',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 77
    def test_frombamfqc(self, ifs):
        ci = [
            "WGBS",
            '-i',
            ifs / 'bam_input' / 'filtered_bam',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 77
    def test_frombamskipqc(self, ifs):
        ci = [
            "WGBS",
            '-i',
            ifs / 'bam_input' / 'filtered_bam',
            '-o',
            'outdir',
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
        assert parseSpOut(_p) == 41

class TestATAC():
    def test_default(self, ifs):
        ci = [
            "ATAC-seq",
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
        assert parseSpOut(_p) == 48
    def test_genrich(self, ifs):
        ci = [
            "ATAC-seq",
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
        assert parseSpOut(_p) == 56
    def test_HMMRATAC(self, ifs):
        ci = [
            "ATAC-seq",
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
        assert parseSpOut(_p) == 55
    def test_sieve(self, ifs):
        ci = [
            "ATAC-seq",
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
        assert parseSpOut(_p) == 48
    def test_frombam(self, ifs):
        ci = [
            "ATAC-seq",
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
        assert parseSpOut(_p) == 84

class TestHIC():
    def test_default(self, ifs):
        ci = [
            "HiC",
            '-i',
            ifs / 'PE',
            '-o',
            'output',
            '--snakemakeOptions',
            SMKOPTS,
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 69
    def test_ice(self, ifs):
        ci = [
            "HiC",
            '-i',
            ifs / 'PE',
            '-o',
            'output',
            '--snakemakeOptions',
            SMKOPTS,
            '--correctionMethod',
            'ICE',
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 75
    def test_trim(self, ifs):
        ci = [
            "HiC",
            '-i',
            ifs / 'PE',
            '-o',
            'output',
            '--snakemakeOptions',
            SMKOPTS,
            '--trim',
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 75
    def test_dpnii(self, ifs):
        ci = [
            "HiC",
            '-i',
            ifs / 'PE',
            '-o',
            'output',
            '--snakemakeOptions',
            SMKOPTS,
            '--enzyme',
            'DpnII',
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 69
    def test_notad(self, ifs):
        ci = [
            "HiC",
            '-i',
            ifs / 'PE',
            '-o',
            'output',
            '--snakemakeOptions',
            SMKOPTS,
            '--noTAD',
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 63
    def test_bwamem2(self, ifs):
        ci = [
            "HiC",
            '-i',
            ifs / 'PE',
            '-o',
            'output',
            '--snakemakeOptions',
            SMKOPTS,
            '--aligner',
            'bwa-mem2',
            ifs / 'org.yaml'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 69

class Testpreprocessing():
    def test_default(self, ifs):
        ci = [
            "preprocessing",
            '-i',
            ifs / 'PE',
            '-o',
            'output',
            '--snakemakeOptions',
            SMKOPTS,
            '--fastqc',
            '--optDedupDist',
            '2500'
        ]
        print(' '.join([str(i) for i in ci]))
        _p = sp.run(ci, capture_output=True, text=True)
        assert _p.returncode == 0
        assert parseSpOut(_p) == 57
    def test_DAG(self, ifs):
        ci = [
            "preprocessing",
            '-i',
            ifs / 'PE',
            '-o',
            'output',
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
        assert parseSpOut(_p) == 57
