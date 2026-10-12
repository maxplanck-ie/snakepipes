"""Download UCSC's RepeatMasker table and convert it to Ensembl-style chrom names / 1-based coords."""

import os
import gzip

from common import UCSC_BASE, log, download_or_warn, ensure_dir


def convert_rmsk_to_ensembl(ucsc_rmsk_gz, out_txt):
    """UCSC -> Ensembl: strip the 'chr' prefix, shift 0-based start to 1-based."""
    with gzip.open(ucsc_rmsk_gz, "rt") as fin, open(out_txt, "w") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            # UCSC rmsk columns: 5 genoName, 6 genoStart, 7 genoEnd, ...
            chrom = fields[5]
            start = int(fields[6]) + 1
            end = int(fields[7])

            if chrom.startswith("chr"):
                chrom = chrom[3:]

            fields[5], fields[6], fields[7] = chrom, str(start), str(end)
            fout.write("\t".join(fields) + "\n")


def fetch_rmsk_txt(ucsc_assembly, out_txt):
    ensure_dir(os.path.dirname(out_txt) or ".")
    url = f"{UCSC_BASE}/{ucsc_assembly}/database/rmsk.txt.gz"
    tmp_gz = out_txt + ".gz"

    if not download_or_warn(url, tmp_gz):
        return False

    convert_rmsk_to_ensembl(tmp_gz, out_txt)
    os.remove(tmp_gz)
    log(f"Ensembl-format RepeatMasker table written: {out_txt}")
    return True
