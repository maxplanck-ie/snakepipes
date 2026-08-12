"""Extract tDNA.fa/rDNA.fa from structural_RNA.bed via bedtools getfasta."""

import os
import re
from pathlib import Path

from common import log, run, which_or_warn


def filter_bed_by_prefix(source_bed, prefix, out_bed):
    """Keep only the BED lines whose column-4 name starts with prefix."""
    selected = 0
    with open(source_bed) as fin, open(out_bed, "w") as fout:
        for line in fin:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            if fields[3].startswith(prefix):
                fout.write(line)
                selected += 1
    log(f"{out_bed}: {selected} interval(s) matching '{prefix}'")
    return selected


def fix_trna_headers(in_fa, out_fa):
    """bedtools -name emits '>name::chr:start-end(strand)'; TEsmall wants '>name:chr:start-end:strand'."""
    fixed_lines = []
    with open(in_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                hdr = line[1:].strip()
                hdr = hdr.replace("::", ":")
                hdr = re.sub(r"\(([-+])\)$", r":\1", hdr)
                hdr = re.sub(r"^sncRNA:tRNA:", "", hdr)
                fixed_lines.append(">" + hdr + "\n")
            else:
                fixed_lines.append(line)
    with open(out_fa, "w") as out:
        out.writelines(fixed_lines)


def extract_named_sequences_from_bed(genome_fa, source_bed, pattern, out_fa, rename_trna=False):
    """Extract intervals whose BED column-4 name matches pattern."""
    tmp_bed = out_fa + ".tmp.bed"
    selected = 0

    with open(source_bed) as fin, open(tmp_bed, "w") as fout:
        for line in fin:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            if pattern(fields[3]):
                fout.write(line)
                selected += 1

    if selected == 0:
        Path(out_fa).write_text("")
        os.remove(tmp_bed)
        log(f"No matching intervals for {out_fa}; wrote empty FASTA")
        return 0

    run(["bedtools", "getfasta", "-s", "-name", "-fi", genome_fa, "-bed", tmp_bed, "-fo", out_fa])

    if rename_trna:
        fixed_lines = []
        with open(out_fa) as fh:
            for line in fh:
                if line.startswith(">"):
                    hdr = line[1:].strip()
                    hdr = hdr.replace("::", ":")
                    hdr = re.sub(r"\(([-+])\)$", r":\1", hdr)
                    hdr = re.sub(r"^sncRNA:tRNA:", "", hdr)
                    fixed_lines.append(">" + hdr + "\n")
                else:
                    fixed_lines.append(line)
        with open(out_fa, "w") as out:
            out.writelines(fixed_lines)

    os.remove(tmp_bed)
    log(f"{out_fa} extracted ({selected} intervals)")
    return selected


def faidx_if_possible(fa_path):
    if not which_or_warn("samtools"):
        Path(fa_path + ".fai").write_text("")
        return False
    if os.path.getsize(fa_path) == 0:
        Path(fa_path + ".fai").write_text("")
        log(f"Empty FASTA; wrote empty index placeholder: {fa_path}.fai")
        return False
    run(["samtools", "faidx", fa_path])
    return True
