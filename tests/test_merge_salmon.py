"""Unit tests for bin/merge_salmon.py (stdlib unittest, no Nextflow needed).

Run:  python -m unittest discover -s tests -v
"""
import importlib.util
import json
import os
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("merge_salmon", os.path.join(HERE, "..", "bin", "merge_salmon.py"))
ms = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ms)

HEADER = "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"


def fake_sample(root, name, tx_rows, gene_rows, pct):
    d = os.path.join(root, name)
    os.makedirs(os.path.join(d, "aux_info"))
    with open(os.path.join(d, "quant.sf"), "w") as fh:
        fh.write(HEADER + "".join(f"{n}\t1000\t850\t{tpm}\t{reads}\n" for n, tpm, reads in tx_rows))
    with open(os.path.join(d, "quant.genes.sf"), "w") as fh:
        fh.write(HEADER + "".join(f"{n}\t1000\t850\t{tpm}\t{reads}\n" for n, tpm, reads in gene_rows))
    with open(os.path.join(d, "aux_info", "meta_info.json"), "w") as fh:
        json.dump({"num_processed": 1000, "percent_mapped": pct}, fh)
    return d


def read_tsv(path):
    with open(path) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh]
    return rows[0], {r[0]: r[1:] for r in rows[1:]}


class TestMerge(unittest.TestCase):
    def test_matrices_and_missing_features(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = fake_sample(tmp, "sA", [("t1", 10, 5), ("t2", 90, 45)], [("g1", 100, 50)], 91.5)
            # sB lacks t2 entirely -> must be filled with 0, not crash
            b = fake_sample(tmp, "sB", [("t1", 100, 20)], [("g1", 100, 20)], 88.0)
            out = os.path.join(tmp, "out")
            os.makedirs(out)
            ms.main([b, a, "--outdir", out])  # deliberately unsorted input

            header, tx_counts = read_tsv(os.path.join(out, "merged_transcript_counts.tsv"))
            self.assertEqual(header, ["transcript_id", "sA", "sB"])  # samples sorted
            self.assertEqual(tx_counts["t2"], ["45.0000", "0.0000"])

            _, gene_tpm = read_tsv(os.path.join(out, "merged_gene_tpm.tsv"))
            self.assertEqual(gene_tpm["g1"], ["100.0000", "100.0000"])

            _, mapping = read_tsv(os.path.join(out, "mapping_summary.tsv"))
            self.assertEqual(mapping["sA"], ["1000", "91.5"])

    def test_bad_header_exits(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = os.path.join(tmp, "bad")
            os.makedirs(d)
            with open(os.path.join(d, "quant.sf"), "w") as fh:
                fh.write("foo\tbar\n")
            with self.assertRaises(SystemExit):
                ms.main([d, "--outdir", tmp])


if __name__ == "__main__":
    unittest.main()
