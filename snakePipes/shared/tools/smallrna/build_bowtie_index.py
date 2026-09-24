"""Build a Bowtie1 index for one FASTA file."""

import os

from common import warn, run, which_or_warn, ensure_dir


def build_bowtie_index(fa_path, out_prefix):
    if not which_or_warn("bowtie-build"):
        warn(f"Skipping bowtie-build for {fa_path}")
        return False
    if not os.path.exists(fa_path) or os.path.getsize(fa_path) == 0:
        warn(f"Skipping bowtie-build for empty/missing FASTA: {fa_path}")
        return False
    ensure_dir(os.path.dirname(out_prefix) or ".")
    run(["bowtie-build", "--threads", "20", fa_path, out_prefix])
    return True
