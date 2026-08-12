"""Download miRBase/piRNAdb annotation, convert to BED6, liftOver if the source build doesn't match."""

import os
import re
import shutil
import zipfile
import tempfile

from common import (
    log, ensure_dir, download_or_raise, liftover,
    normalize_genome, get_species, parse_gff3_attributes, to_ensembl_chrom,
    MIRBASE_URLS, PIRNADB_URLS,
)


def parse_genome_build_header(path, source):
    with open(path) as fh:
        for line in fh:
            if source == "mirbase":
                if line.startswith("# genome-build-id:"):
                    return line.split(":")[-1].strip()
            else:
                if line.startswith("#!genome-build"):
                    return line.split()[-1]
    raise RuntimeError(f"Genome build not found in {path}")


def mirbase_gff_to_bed(gff, feature, outfile):
    """Writes bare Ensembl-style chrom names ('1', not 'chr1') even though
    miRBase's own GFF3 ships 'chr'-prefixed names for human/mouse -- every
    other file in the reference (genome.fa, exon/intron/TE/structural_RNA.bed)
    is Ensembl-style, and there's no separate normalization pass downstream
    for the case where source and target build already match (only the
    liftOver path round-trips through/back from UCSC-style), so this has to
    be correct at the point of writing."""
    with open(gff) as inp, open(outfile, "w") as out:
        for line in inp:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) != 9 or f[2] != feature:
                continue
            attr = parse_gff3_attributes(f[8])
            name = attr.get("Name") or attr.get("ID") or "miRNA"
            chrom = to_ensembl_chrom(f[0])
            out.write(f"{chrom}\t{int(f[3]) - 1}\t{f[4]}\t{name}\t0\t{f[6]}\n")


def pirnadb_gtf_to_bed(gtf, outfile):
    """See mirbase_gff_to_bed: normalizes to bare Ensembl-style chrom names."""
    with open(gtf) as inp, open(outfile, "w") as out:
        for line in inp:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) != 9 or f[2] != "piRNA":
                continue
            m = re.search(r'piRNA_code "([^"]+)"', f[8])
            name = m.group(1) if m else "piRNA"
            chrom = to_ensembl_chrom(f[0])
            out.write(f"{chrom}\t{int(f[3]) - 1}\t{f[4]}\t{name}\t0\t{f[6]}\n")


def download_smallrna_source(source, species, tmp):
    urls = MIRBASE_URLS if source == "mirbase" else PIRNADB_URLS
    ext = "gff3" if source == "mirbase" else "zip"
    outfile = os.path.join(tmp, f"{species}.{source}.{ext}")
    download_or_raise(urls[species], outfile)
    return outfile


def download_smallrna_raw(source, genome, outdir, build_info_name="build_info.txt"):
    """Download + convert to raw (uncollapsed, not-yet-lifted) BED6, and
    record source vs. target build in build_info_name as 'source<TAB>target'.

    build_info_name needs to be distinct per source since mirna and pirna
    both write into annotation/ now -- otherwise one checkpoint's
    build_info.txt would clobber the other's. The actual liftOver call
    happens later, in the Snakefile's checkpoint-gated liftover_bed rule,
    since we can't tell if it's needed until after this download.
    """
    ensure_dir(outdir)
    target_build = normalize_genome(genome)
    species = get_species(genome)
    tmp = tempfile.mkdtemp(prefix=f"{source}_")

    try:
        archive = download_smallrna_source(source, species, tmp)

        if source == "mirbase":
            annotation = archive
            source_build = normalize_genome(parse_genome_build_header(annotation, source))
            mirbase_gff_to_bed(annotation, "miRNA_primary_transcript",
                                os.path.join(outdir, "hairpin.raw.bed"))
            mirbase_gff_to_bed(annotation, "miRNA",
                                os.path.join(outdir, "miRNA.raw.bed"))

        else:  # pirnadb
            with zipfile.ZipFile(archive) as z:
                gtf_name = [x for x in z.namelist() if x.endswith(".gtf")][0]
                z.extract(gtf_name, tmp)
            annotation = os.path.join(tmp, gtf_name)
            source_build = normalize_genome(parse_genome_build_header(annotation, source))
            pirnadb_gtf_to_bed(annotation, os.path.join(outdir, "piRNA_cluster.raw.bed"))

        with open(os.path.join(outdir, build_info_name), "w") as fh:
            fh.write(f"{source_build}\t{target_build}\n")

        if source_build != target_build:
            log(f"{source}: source build {source_build} != target {target_build}; liftOver needed")
        else:
            log(f"{source}: source build matches target ({source_build})")

    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def build_smallrna_annotation(source, genome, outdir):
    """source: 'mirbase' or 'pirnadb'; genome: target build, e.g. 'mm10'.

    Used by the standalone script -- decides on and runs liftOver right
    here instead of going through a checkpoint like the Snakefile does.
    """
    ensure_dir(outdir)
    requested_build = normalize_genome(genome)
    species = get_species(genome)
    tmp = tempfile.mkdtemp(prefix=f"{source}_")

    try:
        archive = download_smallrna_source(source, species, tmp)

        if source == "mirbase":
            annotation = archive
            build = normalize_genome(parse_genome_build_header(annotation, source))

            hairpin = os.path.join(outdir, "hairpin.bed")
            mature = os.path.join(outdir, "miRNA.bed")
            mirbase_gff_to_bed(annotation, "miRNA_primary_transcript", hairpin)
            mirbase_gff_to_bed(annotation, "miRNA", mature)
            files = [hairpin, mature]

        else:  # pirnadb
            with zipfile.ZipFile(archive) as z:
                gtf_name = [x for x in z.namelist() if x.endswith(".gtf")][0]
                z.extract(gtf_name, tmp)
            annotation = os.path.join(tmp, gtf_name)
            build = normalize_genome(parse_genome_build_header(annotation, source))

            pirna = os.path.join(outdir, "piRNA_cluster.bed")
            pirnadb_gtf_to_bed(annotation, pirna)
            files = [pirna]

        if build != requested_build:
            log(f"{source}: lifting {build} -> {requested_build}")
            for f in files:
                lifted = f + ".lifted"
                liftover(f, build, requested_build, lifted)
                os.replace(lifted, f)
        else:
            log(f"{source}: genome already matches ({build})")

    finally:
        shutil.rmtree(tmp, ignore_errors=True)
