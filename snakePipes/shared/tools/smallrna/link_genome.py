"""Symlink an existing genome.fa/.fai into the pipeline's sequence/ directory."""

import os

from common import log, err, ensure_dir


def link_genome_fasta(genome_dir, out_dir):
    ensure_dir(out_dir)
    genome_fa = os.path.join(out_dir, "genome.fa")
    genome_fai = os.path.join(out_dir, "genome.fa.fai")

    for filename in ["genome.fa", "genome.fa.fai"]:
        src = os.path.join(genome_dir, filename)
        dst = os.path.join(out_dir, filename)
        if not os.path.exists(src):
            err(f"Missing file: {src}")
            return None, None
        if os.path.exists(dst) or os.path.islink(dst):
            os.remove(dst)
        os.symlink(src, dst)

    log(f"Linked genome FASTA: {genome_fa}")
    log(f"Linked genome index: {genome_fai}")
    return genome_fa, genome_fai
