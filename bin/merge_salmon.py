#!/usr/bin/env python3
"""
Merge per-sample Salmon outputs into count and TPM matrices.

Reads <sample>/quant.sf (transcript level) and, when present,
<sample>/quant.genes.sf (gene level, written when salmon is run with -g).

Writes:
  merged_transcript_counts.tsv   merged_transcript_tpm.tsv
  merged_gene_counts.tsv         merged_gene_tpm.tsv      (if gene files exist)
  mapping_summary.tsv            per-sample mapping rate from aux_info/meta_info.json

Pure standard library on purpose: it runs in a plain python container.
"""
import argparse
import json
import os
import sys


def read_quant(path):
    """Return {name: (NumReads, TPM)} from a salmon quant(.genes).sf file."""
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_name, i_tpm, i_reads = header.index("Name"), header.index("TPM"), header.index("NumReads")
        except ValueError:
            sys.exit(f"ERROR: unexpected header in {path}: {header}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[i_name]] = (float(f[i_reads]), float(f[i_tpm]))
    return out


def write_matrix(path, samples, data, idx, label):
    features = sorted(set().union(*(d.keys() for d in data.values())))
    with open(path, "w") as fh:
        fh.write(label + "\t" + "\t".join(samples) + "\n")
        for feat in features:
            vals = [data[s].get(feat, (0.0, 0.0))[idx] for s in samples]
            fh.write(feat + "\t" + "\t".join(f"{v:.4f}" for v in vals) + "\n")


def mapping_rate(qdir):
    meta = os.path.join(qdir, "aux_info", "meta_info.json")
    if not os.path.exists(meta):
        return None, None
    with open(meta) as fh:
        m = json.load(fh)
    return m.get("num_processed"), m.get("percent_mapped")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("quant_dirs", nargs="+", help="salmon output directories (dir name = sample name)")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args(argv)

    samples = sorted(os.path.basename(os.path.normpath(d)) for d in a.quant_dirs)
    by_name = {os.path.basename(os.path.normpath(d)): d for d in a.quant_dirs}
    if len(by_name) != len(a.quant_dirs):
        sys.exit("ERROR: duplicate sample directory names")

    tx = {s: read_quant(os.path.join(by_name[s], "quant.sf")) for s in samples}
    write_matrix(os.path.join(a.outdir, "merged_transcript_counts.tsv"), samples, tx, 0, "transcript_id")
    write_matrix(os.path.join(a.outdir, "merged_transcript_tpm.tsv"), samples, tx, 1, "transcript_id")

    gene_files = {s: os.path.join(by_name[s], "quant.genes.sf") for s in samples}
    if all(os.path.exists(p) for p in gene_files.values()):
        genes = {s: read_quant(p) for s, p in gene_files.items()}
        write_matrix(os.path.join(a.outdir, "merged_gene_counts.tsv"), samples, genes, 0, "gene_id")
        write_matrix(os.path.join(a.outdir, "merged_gene_tpm.tsv"), samples, genes, 1, "gene_id")

    with open(os.path.join(a.outdir, "mapping_summary.tsv"), "w") as fh:
        fh.write("sample\treads_processed\tpercent_mapped\n")
        for s in samples:
            n, pct = mapping_rate(by_name[s])
            fh.write(f"{s}\t{n if n is not None else 'NA'}\t{pct if pct is not None else 'NA'}\n")


if __name__ == "__main__":
    main()
