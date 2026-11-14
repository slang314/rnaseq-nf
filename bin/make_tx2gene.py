#!/usr/bin/env python3
"""
Build a transcript->gene table from a GENCODE transcriptome FASTA.

GENCODE headers look like:
  >ENST00000456328.2|ENSG00000290825.1|-|-|DDX11L2-202|DDX11L2|1657|lncRNA|
Field 1 is the transcript ID, field 2 the gene ID and field 6 the gene name.
With `salmon index --gencode` salmon keeps only field 1, so the table below
matches salmon's transcript names.

Usage: make_tx2gene.py transcripts.fa(.gz) [--names gene_names.tsv] > tx2gene.tsv
  --names  also write gene_id<TAB>gene_name (one row per gene, with a header)
"""
import argparse
import gzip
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta")
    ap.add_argument("--names", help="write gene_id<TAB>gene_name to this file")
    a = ap.parse_args(argv)

    opener = gzip.open if a.fasta.endswith(".gz") else open
    n = 0
    names = {}   # gene_id -> gene_name; dicts keep insertion order
    with opener(a.fasta, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            fields = line[1:].strip().split("|")
            if len(fields) < 2:
                sys.exit(f"ERROR: header is not GENCODE-style (no '|'): {line.strip()}\n"
                         "Supply --tx2gene yourself for non-GENCODE references.")
            sys.stdout.write(f"{fields[0]}\t{fields[1]}\n")
            if len(fields) >= 6 and fields[5]:
                names.setdefault(fields[1], fields[5])
            n += 1
    if n == 0:
        sys.exit("ERROR: no FASTA headers found")

    if a.names:
        with open(a.names, "w") as out:
            out.write("gene_id\tgene_name\n")
            for gene_id, name in names.items():
                out.write(f"{gene_id}\t{name}\n")


if __name__ == "__main__":
    main()
