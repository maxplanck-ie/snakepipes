from pathlib import Path
import subprocess as sp
from ruamel.yaml import YAML
import shutil
import gzip
import pytest

def extract_gz(i, o):
  '''
  extracts i into file o
  '''
  with gzip.open(i, 'rb') as f:
    with open(o, 'wb') as of:
      shutil.copyfileobj(f, of)

def createTestData(fp):
    '''
    Sets up fasta, gtf and organism yaml in a factory (fp)
    '''
    # uncompress fna into tmp factory
    fnagz = Path('tests') / 'data' / 'genomes' / 'genome_chr17.fa.gz'
    faout = fp / 'genome.fa'
    extract_gz(fnagz, faout)
    # uncompress gtf into tmp factory
    gtfgz = Path('tests') / 'data' / 'genomes' / 'genes_chr17.gtf.gz'
    gtfout = fp / 'genes.gtf'
    extract_gz(gtfgz, gtfout)
    # Path to STAR index (created in action)
    STARpath = Path('tests') / 'data' / 'mRNA_STAR'
    STARpath = p.resolve()
    
    orgyaml = {
      "genome_size": 94987271 , #we can also extract genome size from STARindex output
      "genome_fasta": faout.as_posix(),
      "star_index": STARpath.as_posix(),
      "genes_gtf" : gtfout.as_posix(),
      "extended_coding_regions_gtf" : "",
      "blacklist_bed": "",
      "ignoreForNormalization": "" 
    }
    # set up yaml
    yaml = YAML()
    yaml.boolean_representation = ['False', 'True']
    with open(fp / 'org.yaml', 'w') as of:
        yaml.dump(orgyaml, of)

@pytest.fixture(scope='session')
def ifs(tmp_path_factory):
  fp = tmp_path_factory.mktemp("data")
  createTestData(fp)
  return fp

class TestmRNAseq:
    def test_mrna(self, ifs):
      org = ifs / 'org.yaml'
      clusterconfig = Path('tests') / 'data' / 'cluster_config.yaml'
      sp.run(
        [
          'mRNA-seq',
          '-i',
          Path('tests') / 'data' / 'mRNA',
          '-o',
          'test_mrna',
          '--clusterConfig',
          clusterconfig,
          org
        ]
      )
      assert Path('test_mrna/mRNA-seq_snakePipes.done').is_file == True

      
