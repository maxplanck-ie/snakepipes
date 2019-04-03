#!/bin/bash
export PATH="$HOME/miniconda/bin:$PATH"
hash -r

cd docs && make html
