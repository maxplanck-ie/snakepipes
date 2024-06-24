import os
from pathlib import Path
import subprocess as sp
import pytest
from ruamel.yaml import YAML



#  for the functionality of tests/test_mrnaseq.py
# * 		function that creates the STAR index from fasta and gtf (note my previous message)
# * 		Create organism yaml with STAR index you just created
# * 		Create a custom config file which sets conda path to be empty
# * 		Have envs created related to mRNA-seq workflow (createIndices functionality)
# * 		implement the actual test



genomedir = Path(__file__).parents[0]/ 'data' / 'genomes'
datadir = Path(__file__).parents[0] / 'data' / 'mRNA'
gtf = genomedir/'genes_chr17.gtf'
fa = genomedir/'genome_chr17.fa'

def STAR_ix():
  assert gtf.exists(), f"GTF file does not exist at {gtf}"
  assert fa.exists(), f"FASTA file does not exist at {fa}"

  sp.run(['STAR' , '--runThreadN' , '8', '--runMode', 'genomeGenerate' , 
          '--genomeDir' , 'STARdir2' , '--genomeFastaFiles' , str(fa), 
          '--sjdbGTFfile',str(gtf), '--sjdbOverhang', '100', '--genomeSAindexNbases' , '12'], check=True)

  return 
       
STAR_ix()


def createYaml(yaml_file):
  yamlPath  = datadir / yaml_file
  _yaml_ = YAML()

  data = {
     "genome_size": 94987271 , #we can also extract genome size from STARindex output
     "genome_fasta": "data/manke/group/schmidth/tests/data/genome_chr17.fa",
     "star_index": "STARdir",
     "genes_gtf" : "data/manke/group/schmidth/tests/data/genomes/genes_chr17.gtf",
     "extended_coding_regions_gtf" : "",
     "blacklist_bed": "",
     "ignoreForNormalization": "" 
     }
  
  with open(yamlPath, 'w') as file: 
     _yaml_.dump(data, file)

  return yamlPath.as_posix()


def config_conda(_config):
  config_path = Path(__file__).parents[0] / 'data' / 'env' / _config
  config_data = { "conda_path": "" }

  with open(config_path, "w") as config_file:
    _yaml = YAML()
    _yaml.dump(config_data, config_file)
  
  return config_path.as_posix()


def setup_env():
  sp.run(['createIndices', '-o', str(genomedir), 
          '--genomeURL', str(fa),  
          '--gtfURL', str(gtf), 
          '--userYAML',  
          'mm10_M19_chr17'], 
            check = True)
  


def setup_files():
    yaml_file = 'mm10_chr17.yaml'
    STAR_index = STAR_ix()
    yaml_path = createYaml(yaml_file)
    _config = 'conda.yaml'
    config= config_conda(_config)

    return STAR_index, yaml_path, config

def test_mrnaseq(setup_files):
    STAR_index , yaml_path, config_path = setup_files

    assert Path(STAR_index).exists(), "STAR index directory does not exist"
    assert Path(yaml_path).exists(), "Organism YAML file does not exist"
    assert Path(config_path).exists(), "Custom config file does not exist"

    setup_env()


# def test_mrna():
#   sp.run(
#     ['mRNA-seq', '-i', '/data/manke/group/schmidth/snakepipes/tests/data/mRNA', '-o', 'TEST', '/data/manke/group/schmidth/snakepipes/tests/mm10_chr17.yaml']
#   )