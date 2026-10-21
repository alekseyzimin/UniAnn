"""Run with: python3 -m unittest discover -s tests -p test_coding_frame.py."""
import math
from pathlib import Path
import subprocess
import tempfile
import unittest


CONVERTER = Path(__file__).resolve().parents[1] / "scripts/convert_coding_frame.pl"
HEADER = ("chrom\tposition_1based\tstrand\tframe0_probability\t"
          "frame1_probability\tframe2_probability\tnoncoding_probability\n")


class CodingFrameTests(unittest.TestCase):
    def convert(self, sequence, probabilities, strand="+", header=False):
        table = HEADER if header else ""
        table += "".join(
            "\t".join(["chr", str(n + 1), strand, *map(str, values)]) + "\n"
            for n, values in enumerate(probabilities)
        )
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "sequence.fa"
            fasta.write_text(f">chr\n{sequence}\n")
            result = subprocess.run(
                ["perl", str(CONVERTER), str(fasta)], input=table,
                capture_output=True, text=True, check=True,
            )
        rows = [line.split("\t") for line in result.stdout.splitlines()]
        self.assertEqual(len(rows), len(sequence))
        self.assertEqual([int(row[0]) for row in rows], list(range(len(sequence))))
        self.assertEqual("".join(row[-1] for row in rows), sequence)
        return rows

    def test_header_does_not_change_rows(self):
        sequence = "CCATGACCC"
        probabilities = [(0.7, 0.1, 0.15, 0.05)] * len(sequence)
        self.assertEqual(self.convert(sequence, probabilities),
                         self.convert(sequence, probabilities, header=True))

    def test_codon_classes_in_all_phases_and_spliced_sequences(self):
        # As in the pipeline, the FASTA is already oriented for each strand.
        for strand in ("+", "-"):
            for start in range(3):
                for intron_length in (0, 4, 5):
                    sequence = "C" * start + "ATGA" + "C" * intron_length + "AACCC"
                    positions = list(range(start, start + 4))
                    positions += list(range(start + 4 + intron_length, len(sequence)))
                    probabilities = [[0.1, 0.1, 0.1, 0.7] for _ in sequence]
                    for local, n in enumerate(positions):
                        probabilities[n] = [0.01, 0.01, 0.01, 0.01]
                        probabilities[n][local % 3] = 0.97
                    rows = self.convert(sequence, probabilities, strand)
                    for local, n in enumerate(positions):
                        with self.subTest(strand=strand, start=start,
                                          intron=intron_length, position=n):
                            scores = list(map(float, rows[n][2:5]))
                            self.assertEqual(scores.index(max(scores)), (n - local) % 3)
                    # Check the current T=3 scaling and coding/noncoding formulas.
                    scaled = [1 / (1 + ((1 - p) / (p + 1e-10)) ** (1 / 3))
                              for p in probabilities[start]]
                    self.assertAlmostEqual(float(rows[start][1]),
                                           math.log(math.e * scaled[3] + 1e-10), delta=0.0051)
                    self.assertAlmostEqual(float(rows[start][2 + start]),
                                           math.log(math.e * scaled[0] * (1 - scaled[3]) + 1e-10),
                                           delta=0.0051)

    def test_stop_masks_follow_the_absolute_frame(self):
        for strand in ("+", "-"):
            for start in range(3):
                for stop in ("TAA", "TAG", "TGA"):
                    with self.subTest(strand=strand, start=start, stop=stop):
                        sequence = "C" * start + stop + "CCC"
                        rows = self.convert(sequence, [(0.7, 0.1, 0.15, 0.05)] * len(sequence), strand)
                        for n, row in enumerate(rows):
                            scores = list(map(float, row[2:5]))
                            if n == start + 2:
                                self.assertEqual(scores[start], -1000000)
                                self.assertEqual(row[5], "-2.00")
                                self.assertEqual(scores.count(-1000000), 1)
                            else:
                                self.assertNotIn(-1000000, scores)
                                self.assertEqual(row[5], row[1])


if __name__ == "__main__":
    unittest.main()
