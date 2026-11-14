#!/usr/bin/env python3
"""
Build a transcript->gene table from a GENCODE transcriptome FASTA.

GENCODE headers look like:
  >ENST00000456328.2|ENSG00000290825.1|-|-|DDX11L2-202|DDX11L2|1657|lncRNA|
Field 1 is the transcript ID and field 2 the gene ID. With `salmon index --gencode`
salmon keeps only field 1, so the table below matches salmon's transcript names.

Usage: make_tx2gene.py transcripts.fa(.gz) > tx2gene.tsv
"""
import gzip
import sys


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    path = sys.argv[1]
    opener = gzip.open if path.endswith(".gz") else open
    n = 0
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            fields = line[1:].strip().split("|")
            if len(fields) < 2:
                sys.exit(f"ERROR: header is not GENCODE-style (no '|'): {line.strip()}\n"
                         "Supply --tx2gene yourself for non-GENCODE references.")
            sys.stdout.write(f"{fields[0]}\t{fields[1]}\n")
            n += 1
    if n == 0:
        sys.exit("ERROR: no FASTA headers found")


if __name__ == "__main__":
    main()
