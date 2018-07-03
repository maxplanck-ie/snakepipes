#!/usr/bin/env python
from setuptools import setup, find_packages
from glob import glob
import os.path

scripts = ['bin/snakePipes']
for d in glob('snakePipes/workflows/*'):
    scripts.append(os.path.join(d, os.path.split(d)[1]))

setup(
    name='snakePipes',
    version='1.0.0',
    scripts=scripts,
    packages=find_packages(),
    include_package_data=True,
    url='https://github.com/maxplanck-ie/snakepipes',
    license='GPL v3',
    description='Snakemake workflows and wrappers for NGS data processing from the MPI-IE',
    install_requires=[
        "snakemake"
    ],
    zip_safe=False,
)
