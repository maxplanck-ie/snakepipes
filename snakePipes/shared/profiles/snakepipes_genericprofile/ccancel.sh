#!/usr/bin/env bash
set -eu
module load slurm
scancel "$@"
