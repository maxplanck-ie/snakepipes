"""Build exon.bed/intron.bed from a GTF, and structural_RNA.bed (rRNA/tRNA/snRNA/scRNA/srpRNA) from the GTF + RepeatMasker table."""

import os
import shutil
from collections import defaultdict

from common import (
    log, write_bed, collapse_bed,
    parse_gtf_attributes, sanitize_repeat_name,
    STRUCTURAL_RNA_FROM_GTF, STRUCTURAL_RNA_FROM_RMSK, TE_CLASSES,
)


# --- exon / intron ---

def read_exons(gtf):
    """Read exon entries -> (flat BED6 records, transcript_id -> exon list)."""
    exons_by_tx = defaultdict(list)
    exon_records = []

    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) != 9 or f[2].lower() != "exon":
                continue

            chrom, start, end, strand, attr = f[0], f[3], f[4], f[6], f[8]
            a = parse_gtf_attributes(attr)

            gene = a.get("gene_id", "NA")
            tx = a.get("transcript_id", "NA")
            start, end = int(start) - 1, int(end)
            exon_no = a.get("exon_number", str(len(exons_by_tx[tx]) + 1))
            name = f"{gene}:{tx}:exon_{exon_no}"

            exon_records.append((chrom, start, end, name, 0, strand))
            exons_by_tx[tx].append((chrom, start, end, strand, gene))

    return exon_records, exons_by_tx


def write_raw_exon_bed(gtf, outfile):
    """Uncollapsed exon BED6; sort+merge happens later, in its own conda env (see Snakefile's collapse_bed)."""
    records, _ = read_exons(gtf)
    write_bed(records, outfile)


def make_exon_bed(gtf, outfile):
    records, _ = read_exons(gtf)
    raw = outfile + ".raw"
    write_bed(records, raw)
    collapse_bed(raw, outfile)
    os.remove(raw)


def write_raw_intron_bed(gtf, outfile):
    """Uncollapsed intron BED6 -- see write_raw_exon_bed."""
    _, exons_by_tx = read_exons(gtf)
    records = []

    for tx, exons in exons_by_tx.items():
        exons.sort(key=lambda x: x[1])
        chrom, _, _, strand, gene = exons[0]

        for i in range(len(exons) - 1):
            start, end = exons[i][2], exons[i + 1][1]
            if start >= end:
                continue
            name = f"{gene}:{tx}:intron_{i + 1}"
            records.append((chrom, start, end, name, 0, strand))

    write_bed(records, outfile)


def make_intron_bed(gtf, outfile):
    _, exons_by_tx = read_exons(gtf)
    records = []

    for tx, exons in exons_by_tx.items():
        exons.sort(key=lambda x: x[1])
        chrom, _, _, strand, gene = exons[0]

        for i in range(len(exons) - 1):
            start, end = exons[i][2], exons[i + 1][1]
            if start >= end:
                continue
            name = f"{gene}:{tx}:intron_{i + 1}"
            records.append((chrom, start, end, name, 0, strand))

    raw = outfile + ".raw"
    write_bed(records, raw)
    collapse_bed(raw, outfile)
    os.remove(raw)


# --- structural RNA ---

def extract_structural_rna_from_gtf(gtf, outfile):
    """rRNA/tRNA/snRNA transcript-level entries from an Ensembl GTF.
    Name: sncRNA:<class>:<gene_id>:<transcript_id>_copy<N>."""
    counts = defaultdict(int)

    with open(outfile, "w") as out, open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2].lower() != "transcript":
                continue

            chrom, start, end, strand = fields[0], fields[3], fields[4], fields[6]
            attrs = parse_gtf_attributes(fields[8])

            gene_id = attrs.get("gene_id", "NA")
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                continue

            biotypes = [
                attrs.get("gene_type"), attrs.get("gene_biotype"),
                attrs.get("transcript_type"), attrs.get("transcript_biotype"),
            ]
            rna_class = None
            for bt in biotypes:
                if bt and bt.lower() in STRUCTURAL_RNA_FROM_GTF:
                    rna_class = STRUCTURAL_RNA_FROM_GTF[bt.lower()]
                    break
            if not rna_class:
                continue

            counts[(rna_class, gene_id)] += 1
            copy_n = counts[(rna_class, gene_id)]
            name = f"sncRNA:{rna_class}:{gene_id}:{transcript_id}_copy{copy_n}"
            out.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")

    log(f"GTF structural RNA entries: {sum(counts.values())}")


def extract_structural_rna_from_rmsk(rmsk_txt, outfile):
    """tRNA/scRNA/srpRNA/rRNA repeat entries from the (Ensembl-coordinate) RepeatMasker table."""
    counts = defaultdict(int)

    with open(outfile, "w") as out, open(rmsk_txt) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip().split()
            if len(fields) < 17:
                continue

            chrom = fields[5]
            start = int(fields[6]) + 1
            end = int(fields[7])
            strand = fields[9]
            rep_name = sanitize_repeat_name(fields[10])
            rep_class = sanitize_repeat_name(fields[11]).lower()

            if strand == "C":
                strand = "-"
            elif strand not in {"+", "-"}:
                strand = "."

            if rep_class not in STRUCTURAL_RNA_FROM_RMSK:
                continue
            rna_type = STRUCTURAL_RNA_FROM_RMSK[rep_class]

            counts[(rna_type, rep_name)] += 1
            copy_n = counts[(rna_type, rep_name)]
            name = f"sncRNA:{rna_type}:{rep_name}:{rep_name}_copy{copy_n}"
            out.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")

    log(f"RMSK structural RNA entries: {sum(counts.values())}")


def extract_te_from_rmsk(rmsk_txt, outfile, te_classes=TE_CLASSES):
    """TE entries from the (Ensembl-coordinate) RepeatMasker table, restricted to te_classes.

    Name column comes out as class:family:subfamily:instance, e.g.
    LTR:Gypsy:IDEFIX_LTR:IDEFIX_LTR_copy1 -- subfamily is repName, and
    instance is subfamily plus a running _copy<N> counter per
    (class, family, subfamily) so each genomic copy gets a unique name.

    Not collapsed like exon/intron/structural_RNA -- every row here is
    already a distinct RepeatMasker call, so there's nothing to merge.
    """
    wanted = {c.lower() for c in te_classes}
    counts = defaultdict(int)

    with open(outfile, "w") as out, open(rmsk_txt) as fh:
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip().split()
            if len(fields) < 17:
                continue

            chrom = fields[5]
            start = int(fields[6])
            end = int(fields[7])
            strand = fields[9]
            rep_name = sanitize_repeat_name(fields[10])      # subfamily
            rep_class = sanitize_repeat_name(fields[11])      # class
            rep_family = sanitize_repeat_name(fields[12])     # family

            if strand == "C":
                strand = "-"
            elif strand not in {"+", "-"}:
                strand = "."

            if rep_class.lower() not in wanted:
                continue

            counts[(rep_class, rep_family, rep_name)] += 1
            copy_n = counts[(rep_class, rep_family, rep_name)]
            instance = f"{rep_name}_copy{copy_n}"
            name = f"{rep_class}:{rep_family}:{rep_name}:{instance}"

            out.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")

    log(f"TE entries ({','.join(sorted(te_classes))}): {sum(counts.values())}")


def write_raw_structural_rna(gtf, rmsk, outfile):
    """GTF + RepeatMasker structural RNA entries, concatenated but not yet
    collapsed -- collapsing happens later, in its own conda env."""
    tmp_gtf = outfile + ".gtf.tmp"
    tmp_rmsk = outfile + ".rmsk.tmp"

    extract_structural_rna_from_gtf(gtf, tmp_gtf)
    extract_structural_rna_from_rmsk(rmsk, tmp_rmsk)

    with open(outfile, "w") as out:
        for f in [tmp_gtf, tmp_rmsk]:
            with open(f) as fh:
                shutil.copyfileobj(fh, out)

    os.remove(tmp_gtf)
    os.remove(tmp_rmsk)


def build_structural_rna(gtf, rmsk, outfile):
    tmp_gtf = outfile + ".gtf.tmp"
    tmp_rmsk = outfile + ".rmsk.tmp"

    extract_structural_rna_from_gtf(gtf, tmp_gtf)
    extract_structural_rna_from_rmsk(rmsk, tmp_rmsk)

    with open(outfile, "w") as out:
        for f in [tmp_gtf, tmp_rmsk]:
            with open(f) as fh:
                shutil.copyfileobj(fh, out)

    os.remove(tmp_gtf)
    os.remove(tmp_rmsk)
    collapse_bed(outfile, outfile)
    log(f"Created structural RNA BED: {outfile}")
