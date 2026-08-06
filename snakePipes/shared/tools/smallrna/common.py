"""Shared constants and helpers used by the other stage modules. Not runnable on its own."""

import os
import re
import shutil
import subprocess
import tempfile
import urllib.request
import urllib.error
from pathlib import Path


# --- constants ---

UCSC_BASE = "https://hgdownload.soe.ucsc.edu/goldenPath"

MIRBASE_URLS = {
    "human": "https://www.mirbase.org/download/hsa.gff3",
    "mouse": "https://www.mirbase.org/download/mmu.gff3",
    "drosophila": "https://www.mirbase.org/download/dme.gff3",
}

PIRNADB_URLS = {
    "human": "https://www.pirnadb.org/download/downloadarchive/gff_gtf/pirnadb.v1_7_6.hg38.gtf.zip",
    "mouse": "https://www.pirnadb.org/download/downloadarchive/gff_gtf/pirnadb.v1_7_6.mm10.gtf.zip",
    "drosophila": "https://www.pirnadb.org/download/downloadarchive/gff_gtf/pirnadb.v1_7_6.dm6.gtf.zip",
}

# alias -> canonical UCSC assembly name (also used directly for rmsk/chain downloads)
GENOME_ALIASES = {
    "hg19": "hg19", "GRCh37": "hg19",
    "hg38": "hg38", "GRCh38": "hg38",
    "mm9": "mm9", "NCBI37": "mm9",
    "mm10": "mm10", "GRCm38": "mm10",
    "mm39": "mm39", "GRCm39": "mm39",
    "dm3": "dm3", "Release_5": "dm3",
    "dm6": "dm6", "Release_6": "dm6",
}

SPECIES_BY_GENOME = {
    "human": ["hg19", "hg38", "GRCh37", "GRCh38"],
    "mouse": ["mm9", "mm10", "mm39", "NCBI37", "GRCm38", "GRCm39"],
    "drosophila": ["dm3", "dm6", "Release_5", "Release_6"],
}

STRUCTURAL_RNA_FROM_GTF = {
    "rrna": "rRNA",
    "rrna_pseudogene": "rRNA",
    "trna": "tRNA",
    "snrna": "snRNA",
}

STRUCTURAL_RNA_FROM_RMSK = {
    "trna": "tRNA",
    "scrna": "scRNA",
    "srprna": "srpRNA",
    "rrna": "rRNA",
}

# repClass values that count as TEs for TE.bed -- Simple_repeat/Low_complexity/
# etc are left out on purpose
TE_CLASSES = ["DNA", "LINE", "LTR", "RC", "Retroposon", "RNA", "Satellite", "SINE", "Unknown"]


# --- logging + subprocess ---

def log(msg):
    print(f"[INFO] {msg}", flush=True)


def warn(msg):
    print(f"[WARN] {msg}", flush=True)


def err(msg):
    import sys
    print(f"[ERROR] {msg}", file=sys.stderr, flush=True)


def run(cmd, check=True):
    log("RUN: " + " ".join(cmd))
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if check and result.returncode != 0:
        err(result.stderr.strip())
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return result


def which_or_warn(tool):
    path = shutil.which(tool)
    if not path:
        warn(f"{tool} not found in PATH")
    return path


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def write_empty_bed6(path):
    Path(path).write_text("")


# --- downloads ---

def download_or_raise(url, outfile):
    """Raises on failure -- for files the pipeline can't proceed without (miRBase/piRNAdb, chains)."""
    log(f"Downloading {url}")
    try:
        urllib.request.urlretrieve(url, outfile)
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"Download failed ({e.code}): {url}")


def download_or_warn(url, dest):
    """Returns True/False instead of raising -- for stuff we can limp along without (RepeatMasker)."""
    log(f"Downloading: {url}")
    try:
        with urllib.request.urlopen(url) as response, open(dest, "wb") as out:
            shutil.copyfileobj(response, out)
        return True
    except Exception as e:
        warn(f"Failed download: {url} ({e})")
        return False


# --- genome build / liftover ---

def normalize_genome(genome):
    if genome not in GENOME_ALIASES:
        raise RuntimeError(f"Unsupported genome: {genome}")
    return GENOME_ALIASES[genome]


def get_species(genome):
    for species, genomes in SPECIES_BY_GENOME.items():
        if genome in genomes:
            return species
    raise RuntimeError(f"Cannot determine species from {genome}")


def download_chain(source, target, tmp):
    chain = f"{source}To{target[0].upper() + target[1:]}.over.chain.gz"
    url = f"{UCSC_BASE}/{source}/liftOver/{chain}"
    outfile = os.path.join(tmp, chain)
    try:
        download_or_raise(url, outfile)
    except Exception:
        raise RuntimeError(
            f"No UCSC chain available: {source} -> {target}. "
            "Cannot automatically convert coordinates."
        )
    return outfile


def _bed_uses_chr_prefix(bed_path):
    # peek at the first data line -- 'chr1' (UCSC-style) or '1' (Ensembl-style)?
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            return line.split("\t", 1)[0].startswith("chr")
    return True  # empty file, doesn't matter either way


def _rewrite_chrom_column(src, dst, convert):
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            fields[0] = convert(fields[0])
            fout.write("\t".join(fields) + "\n")


def _to_ucsc_chrom(chrom):
    if chrom.startswith("chr"):
        return chrom
    if chrom in ("MT", "mt", "M"):
        return "chrM"
    return "chr" + chrom


def _to_ensembl_chrom(chrom):
    if not chrom.startswith("chr"):
        return chrom
    chrom = chrom[3:]
    return "MT" if chrom == "M" else chrom


def liftover(infile, source, target, outfile):
    """LiftOver a BED6 file between UCSC assemblies.

    UCSC chain files want 'chr1'-style names. miRBase/piRNAdb ship
    Ensembl-style names ('1', 'MT') with no prefix, and feeding those to
    liftOver as-is just drops almost everything into .unmapped. So we
    check which style the input uses, add 'chr' for the liftOver call,
    then strip it back off so the output matches what came in.
    """
    tmp = tempfile.mkdtemp(prefix="chain_")
    try:
        chain = download_chain(source, target, tmp)

        input_has_chr_prefix = _bed_uses_chr_prefix(infile)
        lift_input = infile
        if not input_has_chr_prefix:
            lift_input = os.path.join(tmp, "input.chr.bed")
            _rewrite_chrom_column(infile, lift_input, _to_ucsc_chrom)

        lift_output = outfile if input_has_chr_prefix else os.path.join(tmp, "lifted.chr.bed")
        unmapped = lift_output + ".unmapped"

        run(["liftOver", lift_input, chain, lift_output, unmapped])

        if os.path.exists(unmapped):
            n_unmapped = sum(1 for l in open(unmapped) if l.strip() and not l.startswith("#"))
            if n_unmapped:
                warn(f"liftOver ({source} -> {target}): {n_unmapped} interval(s) could not be "
                     f"mapped and were dropped from {os.path.basename(outfile)}")
            os.remove(unmapped)

        if not input_has_chr_prefix:
            _rewrite_chrom_column(lift_output, outfile, _to_ensembl_chrom)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --- gtf/gff3 attribute parsing ---

def parse_gff3_attributes(text):
    """key=value;key=value style, used by miRBase GFF3."""
    out = {}
    for item in text.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out


def parse_gtf_attributes(attr):
    """key "value"; style, used by Ensembl GTFs."""
    attrs = {}
    for item in attr.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        match = re.match(r'(\S+)\s+(.+)', item)
        if match:
            attrs[match.group(1)] = match.group(2).strip().strip('"')
    return attrs


def sanitize_repeat_name(name):
    name = name.strip()
    name = re.sub(r"\s+", "_", name)
    name = re.sub(r"[^A-Za-z0-9_.:+-]", "_", name)
    return name


# --- bed helpers ---

def write_bed(records, outfile):
    with open(outfile, "w") as out:
        for r in records:
            out.write("\t".join(map(str, r)) + "\n")


def collapse_bed(raw, output):
    """Sort + merge identical intervals (same chrom/start/end/strand), concatenating their names."""
    sorted_bed = raw + ".sorted"
    grouped_bed = raw + ".grouped"

    subprocess.run(["bedtools", "sort", "-i", raw], stdout=open(sorted_bed, "w"), check=True)
    subprocess.run(
        ["bedtools", "groupby", "-i", sorted_bed, "-g", "1,2,3,6", "-c", "4", "-o", "collapse"],
        stdout=open(grouped_bed, "w"), check=True,
    )

    with open(grouped_bed) as fin, open(output, "w") as fout:
        for line in fin:
            chrom, start, end, strand, names = line.rstrip().split("\t")
            fout.write(f"{chrom}\t{start}\t{end}\t{names}\t0\t{strand}\n")

    os.remove(sorted_bed)
    os.remove(grouped_bed)
