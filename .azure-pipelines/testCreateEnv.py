#!/usr/bin/env python

import subprocess as sp
import sys
sys.path.append("..")
from snakePipes import common_functions as cf
keys = []
index = len(cf.set_env_yamls().items())/10
index = round(index)
for i in range(index):
    keys.append(list(cf.set_env_yamls())[i])
envs_to_test = " ".join(keys)
sp.check_output( "snakePipes createEnvs --only {} --force".format(envs_to_test), shell =True)
