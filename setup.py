#!/usr/bin/env python
from setuptools import setup, find_packages
from glob import glob
import os.path

# Set __version__
exec(open('snakePipes/__init__.py').read())

scripts = ['bin/snakePipes']
for d in glob('snakePipes/workflows/*'):
    scripts.append(os.path.join(d, os.path.split(d)[1]))

setup(
    name='snakePipes',
    version=__version__,  # noqa: F821
    scripts=scripts,
    packages=find_packages(),
    include_package_data=True,
    url='https://github.com/maxplanck-ie/snakepipes',
    license='GPL v3',
    description='Snakemake workflows and wrappers for NGS data processing from the MPI-IE',
    install_requires=[
        "snakemake >= 5.2.0",
        "psutil",
        "pandas",
        "fuzzywuzzy"
    ],
    zip_safe=False,
)
