import os
import subprocess
import sys
import tempfile
import unittest

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin", "make_tx2gene.py")


class TestTx2Gene(unittest.TestCase):
    def run_on(self, text):
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
            fh.write(text)
        try:
            return subprocess.run([sys.executable, SCRIPT, fh.name], capture_output=True, text=True)
        finally:
            os.unlink(fh.name)

    def test_gencode_header(self):
        r = self.run_on(">ENST00000456328.2|ENSG00000290825.1|-|-|DDX11L2-202|DDX11L2|1657|lncRNA|\nACGT\n")
        self.assertEqual(r.returncode, 0)
        self.assertEqual(r.stdout, "ENST00000456328.2\tENSG00000290825.1\n")

    def test_gene_names_file(self):
        fa = (">ENST1.1|ENSG1.1|-|-|DUSP1-201|DUSP1|100|protein_coding|\nACGT\n"
              ">ENST2.1|ENSG1.1|-|-|DUSP1-202|DUSP1|90|protein_coding|\nACGT\n"
              ">ENST3.1|ENSG2.1|-|-|KLF15-201|KLF15|80|protein_coding|\nACGT\n")
        with tempfile.TemporaryDirectory() as d:
            fasta, names = os.path.join(d, "tx.fa"), os.path.join(d, "names.tsv")
            with open(fasta, "w") as fh:
                fh.write(fa)
            r = subprocess.run([sys.executable, SCRIPT, fasta, "--names", names], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(len(r.stdout.splitlines()), 3)   # tx2gene: one row per transcript
            with open(names) as fh:                           # names: one row per gene
                self.assertEqual(fh.read(), "gene_id\tgene_name\nENSG1.1\tDUSP1\nENSG2.1\tKLF15\n")

    def test_non_gencode_header_fails_clearly(self):
        r = self.run_on(">tx1\nACGT\n")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not GENCODE-style", r.stderr)


if __name__ == "__main__":
    unittest.main()
