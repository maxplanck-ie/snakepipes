#!/bin/bash
export PATH="$HOME/miniconda/bin:$PATH"
hash -r
python -m pip install --no-deps --ignore-installed .

cd docs && make html
