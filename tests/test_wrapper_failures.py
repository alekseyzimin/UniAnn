"""Run wrapper regressions with Bash, Perl and small command stubs."""

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class WrapperTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="uniann wrapper ")
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        shutil.copy2(ROOT / "scripts" / "uniann.sh", self.directory / "uniann.sh")
        self.stub("convert_coding_frame.pl", 'cat "$1" -\n')
        self.stub("uniann", "touch native-called\ncat out.ps.txt\n")
        self.stub("gffread", "cat\n")
        self.fasta = self.directory / "genome input.fa"
        self.coding = self.directory / "coding probabilities.tsv"
        self.scores = self.directory / "site scores.tsv"
        self.emissions = self.directory / "out.ps.txt"
        self.output = self.directory / (self.fasta.name + ".uniann.gff")
        self.fasta.write_text(">chr\nCCC\n")
        self.coding.write_text("current coding probabilities\n")
        self.scores.write_text("".join(
            f"chr\t2\t+\t{site}\tNN\t0.9\n"
            for site in ("donor", "acceptor", "start", "stop")
        ))

    def stub(self, name, body):
        executable = self.directory / name
        executable.write_text("#!/bin/bash\n" + body)
        executable.chmod(0o755)

    def run_wrapper(self):
        environment = dict(os.environ, PATH=str(self.directory) + os.pathsep + os.environ["PATH"])
        return subprocess.run(
            ["bash", str(self.directory / "uniann.sh"), "-f", self.fasta.name,
             "-c", self.coding.name, "-s", self.scores.name],
            cwd=self.directory, env=environment, capture_output=True, text=True,
        )

    def assert_failed(self, result, before_decoder=False):
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.output.exists())
        self.assertNotIn("Output gff file is", result.stdout)
        if before_decoder:
            self.assertFalse((self.directory / "native-called").exists())

    def test_success_and_rerun_use_current_inputs(self):
        self.emissions.write_text("stale emissions\n")
        for sequence, coding in (("CCC", "first"), ("AAA", "updated")):
            self.fasta.write_text(f">chr\n{sequence}\n")
            self.coding.write_text(coding + "\n")
            result = self.run_wrapper()
            self.assertEqual(result.returncode, 0, result.stderr)
            expected = self.fasta.read_text() + self.coding.read_text()
            self.assertEqual(self.emissions.read_text(), expected)
            self.assertEqual(self.output.read_text(), expected)
            self.assertIn("Output gff file is", result.stdout)
            for site_file in ("out.gt.txt", "out.ag.txt", "out.atg.txt", "out.stop.txt"):
                coordinate, score = (self.directory / site_file).read_text().split()
                self.assertEqual(coordinate, "1")
                self.assertAlmostEqual(float(score), 1000, delta=1)
            self.assertFalse((self.directory / "out.ps.txt.tmp").exists())

    def test_conversion_failure_preserves_emissions_and_stops_decoder(self):
        self.emissions.write_text("old emissions\n")
        self.stub("convert_coding_frame.pl", "echo partial-emissions\nexit 23\n")
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 23, result.stderr)
        self.assertEqual(self.emissions.read_text(), "old emissions\n")
        self.assert_failed(result, before_decoder=True)

    def test_invalid_scores_stop_decoder(self):
        self.scores.write_text("chr\t2\t+\tdonor\tGT\t0\n")
        result = self.run_wrapper()
        self.assertIn("Incorrect scores", result.stderr)
        self.assert_failed(result, before_decoder=True)

    def test_each_site_output_failure_stops_decoder(self):
        for name in ("out.gt.txt", "out.ag.txt", "out.atg.txt", "out.stop.txt"):
            with self.subTest(name=name):
                path = self.directory / name
                path.unlink(missing_ok=True)
                path.mkdir()
                result = self.run_wrapper()
                self.assert_failed(result, before_decoder=True)
                path.rmdir()

    def test_native_failure_is_not_published(self):
        self.stub("uniann", "touch native-called\necho native-failure >&2\nexit 17\n")
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 17, result.stderr)
        self.assertIn("native-failure", (self.directory / "out.err").read_text())
        self.assert_failed(result)

    def test_either_gffread_failure_is_not_published(self):
        for stage in ("first", "last"):
            with self.subTest(stage=stage):
                condition = '[[ "$1" == --tlf ]]' if stage == "first" else '[[ "$1" != --tlf ]]'
                self.stub("gffread", f"cat\nif {condition}; then exit 19; fi\n")
                result = self.run_wrapper()
                self.assertEqual(result.returncode, 19, result.stderr)
                self.assert_failed(result)

    def test_emission_rename_failure_stops_decoder(self):
        self.stub("mv", "exit 21\n")
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 21, result.stderr)
        self.assert_failed(result, before_decoder=True)

    def test_final_rename_failure_is_not_announced(self):
        self.stub("mv", 'if [[ "$1" == out.ps.txt.tmp ]]; then /bin/mv "$@"; else exit 22; fi\n')
        result = self.run_wrapper()
        self.assertEqual(result.returncode, 22, result.stderr)
        self.assert_failed(result)


if __name__ == "__main__":
    unittest.main()
