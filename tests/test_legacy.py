"""Protect the default decoder's GFF and complete DP/BT dump byte for byte.

Goldens were captured from unmodified upstream commit 0d89ff0. The GFF is
stored as readable text; the much larger diagnostic dump is stored as SHA256.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from fixtures import cases, write_inputs

ROOT = Path(__file__).resolve().parents[1]


class LegacyRegressionTests(unittest.TestCase):
    def test_upstream_output(self):
        expected = json.loads((ROOT / 'tests' / 'legacy_expected.json').read_text())
        for name, sequence, random_scores in cases():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                files = write_inputs(directory, sequence, random_scores)
                result = subprocess.run([str(ROOT / 'src' / 'uniann'), *map(str, files)],
                                        capture_output=True, check=True)
                self.assertEqual(result.stdout.decode(), expected[name]['gff'])
                self.assertEqual(hashlib.sha256(result.stderr).hexdigest(),
                                 expected[name]['dump_sha256'])
                quiet = subprocess.run([str(ROOT / 'src' / 'uniann'), *map(str, files), '--no-dp-dump'],
                                       capture_output=True, check=True)
                self.assertEqual(quiet.stdout, result.stdout)
                self.assertEqual(quiet.stderr, b'')

    def test_rejected_inputs(self):
        missing = subprocess.run([str(ROOT / 'src' / 'uniann')], capture_output=True)
        self.assertEqual(missing.returncode, 1)
        self.assertIn(b'atg.txt stop.txt', missing.stderr)
        for filename, bad_line in [('emissions.txt', '-1 0 0 0 0 0\n'),
                                   ('emissions.txt', '0 0 1\n'), ('gt.txt', '9 1\n')]:
            with self.subTest(filename=filename, line=bad_line), tempfile.TemporaryDirectory() as directory:
                files = write_inputs(directory, 'AAA')
                (Path(directory) / filename).write_text(bad_line)
                result = subprocess.run([str(ROOT / 'src' / 'uniann'), *map(str, files)], capture_output=True)
                self.assertEqual(result.returncode, 1)
                self.assertIn(b'Invalid', result.stderr)
            unknown = subprocess.run([str(ROOT / 'src' / 'uniann'), *map(str, files), '--unknown'], capture_output=True)
            self.assertEqual(unknown.returncode, 1)
            self.assertIn(b'Unknown option', unknown.stderr)


if __name__ == '__main__':
    unittest.main()
