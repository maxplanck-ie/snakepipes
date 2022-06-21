#!/usr/bin/env python

import subprocess as sp
import sys
import argparse

sys.path.append("..")
from snakePipes import common_functions as cf


def parse_args():
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('Required Arguments')
    required.add_argument("-d",
                          dest="envs_chunk",
                          help="chunk of env to be used when running `snakePipes createEnvs`",
                          required=True)
    return parser


parser = parse_args()
args = parser.parse_args()

keys = []
length = len(cf.set_env_yamls().items())
index = length / 3
index = round(index)
if args.envs_chunk == str(1):
    envs = range(0, index)
elif args.envs_chunk == str(2):
    envs = range(index, 2 * index)
elif args.envs_chunk == str(3):
    envs = range(2 * index, length)

for i in range(envs):
    keys.append(list(cf.set_env_yamls())[i])
envs_to_test = " ".join(keys)
sp.check_output("snakePipes createEnvs --only {} --force".format(envs_to_test), shell=True)
