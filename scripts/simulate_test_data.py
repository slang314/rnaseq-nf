#!/usr/bin/env python3
"""
Generate a tiny, fully synthetic bulk RNA-seq dataset with KNOWN answers.

Why synthetic? Because you know the ground truth, so the pipeline test can
check that the biology came out right (simulated fold changes are recovered),
not just that the pipeline finished.

Outputs (in --outdir):
  transcripts.fa        mini transcriptome (8 genes, ~11 transcripts)
  tx2gene.tsv           transcript -> gene map (2 columns, no header)
  truth.tsv             simulated relative abundance per transcript per condition
  <sample>_R1/R2.fastq.gz  paired-end 2x100 reads, 4 samples (2 ctrl, 2 trt)
  samplesheet.csv       pipeline input (paths relative to the samplesheet)

Ground truth baked in:
  GENE_UP    is 4x higher in "treated"
  GENE_DOWN  is 4x lower  in "treated"
  all other genes unchanged

Usage:
  python3 scripts/simulate_test_data.py --outdir assets/test_data
"""
import argparse
import gzip
import os
import random

ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"  # Illumina TruSeq read-1 adapter
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(rng, seq, err):
    out = []
    for b in seq:
        if rng.random() < err:
            out.append(rng.choice([x for x in "ACGT" if x != b]))
        else:
            out.append(b)
    return "".join(out)


def build_transcriptome(rng):
    """8 genes; some with alternative isoforms that share exons (realistic
    multi-mapping that salmon has to resolve)."""
    genes = ["GENE_UP", "GENE_DOWN", "GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E", "GENE_F"]
    tx = {}       # tx_id -> seq
    tx2gene = {}  # tx_id -> gene
    for g in genes:
        exons = [rand_seq(rng, rng.randint(250, 600)) for _ in range(rng.randint(3, 5))]
        n_iso = 1 if g in ("GENE_UP", "GENE_DOWN") else rng.choice([1, 2, 2])
        for i in range(n_iso):
            if i == 0:
                parts = exons
            else:  # skip one internal exon -> alternative isoform
                skip = rng.randrange(1, len(exons) - 1)
                parts = [e for j, e in enumerate(exons) if j != skip]
            tid = f"{g}-T{i + 1}"
            tx[tid] = "".join(parts)
            tx2gene[tid] = g
    return tx, tx2gene


def abundances(rng, tx, tx2gene):
    base = {t: rng.uniform(1, 10) for t in tx}
    ctrl, trt = dict(base), dict(base)
    for t, g in tx2gene.items():
        if g == "GENE_UP":
            trt[t] = base[t] * 4
        elif g == "GENE_DOWN":
            trt[t] = base[t] / 4
    return {"control": ctrl, "treated": trt}


def simulate_sample(rng, tx, abund, n_pairs, read_len=100, err=0.005):
    # Read sampling probability ~ abundance * length (like real RNA-seq)
    tids = list(tx)
    weights = [abund[t] * len(tx[t]) for t in tids]
    r1s, r2s = [], []
    qual = "I" * read_len
    for i in range(n_pairs):
        t = rng.choices(tids, weights)[0]
        s = tx[t]
        # 3% of fragments are short -> read-through into adapter (fastp trims these)
        frag = rng.randint(60, 95) if rng.random() < 0.03 else max(read_len, int(rng.gauss(250, 30)))
        frag = min(frag, len(s))
        start = rng.randint(0, len(s) - frag)
        f = s[start:start + frag]
        if rng.random() < 0.5:  # unstranded library
            f = revcomp(f)
        r1 = (f + ADAPTER + "A" * read_len)[:read_len]
        r2 = (revcomp(f) + ADAPTER_R2 + "A" * read_len)[:read_len]
        name = f"@sim{i}_{t}"
        r1s.append(f"{name}/1\n{mutate(rng, r1, err)}\n+\n{qual}\n")
        r2s.append(f"{name}/2\n{mutate(rng, r2, err)}\n+\n{qual}\n")
    return r1s, r2s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="assets/test_data")
    ap.add_argument("--pairs", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    rng = random.Random(a.seed)

    tx, tx2gene = build_transcriptome(rng)
    ab = abundances(rng, tx, tx2gene)

    with open(os.path.join(a.outdir, "transcripts.fa"), "w") as fh:
        for t, s in tx.items():
            fh.write(f">{t}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    with open(os.path.join(a.outdir, "tx2gene.tsv"), "w") as fh:
        for t, g in tx2gene.items():
            fh.write(f"{t}\t{g}\n")
    with open(os.path.join(a.outdir, "truth.tsv"), "w") as fh:
        fh.write("transcript_id\tgene_id\tcontrol\ttreated\n")
        for t in tx:
            fh.write(f"{t}\t{tx2gene[t]}\t{ab['control'][t]:.3f}\t{ab['treated'][t]:.3f}\n")

    samples = [("ctrl_rep1", "control"), ("ctrl_rep2", "control"),
               ("trt_rep1", "treated"), ("trt_rep2", "treated")]
    rows = ["sample,fastq_1,fastq_2,condition"]
    for name, cond in samples:
        r1, r2 = simulate_sample(rng, tx, ab[cond], a.pairs)
        for reads, tag in ((r1, "R1"), (r2, "R2")):
            # mtime=0 -> byte-identical output on every run (reproducible test data)
            with gzip.GzipFile(os.path.join(a.outdir, f"{name}_{tag}.fastq.gz"), "wb", mtime=0) as fh:
                fh.write("".join(reads).encode())
        rows.append(f"{name},{name}_R1.fastq.gz,{name}_R2.fastq.gz,{cond}")
    with open(os.path.join(a.outdir, "samplesheet.csv"), "w") as fh:
        fh.write("\n".join(rows) + "\n")
    print(f"Wrote {len(tx)} transcripts, {len(samples)} samples x {a.pairs} pairs to {a.outdir}")


if __name__ == "__main__":
    main()
