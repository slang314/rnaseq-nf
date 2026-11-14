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

    def test_non_gencode_header_fails_clearly(self):
        r = self.run_on(">tx1\nACGT\n")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not GENCODE-style", r.stderr)


if __name__ == "__main__":
    unittest.main()
